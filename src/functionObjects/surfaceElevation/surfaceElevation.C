/*---------------------------------------------------------------------------*\
            Copyright (c) 2026, Asim Onder
-------------------------------------------------------------------------------
License
    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.
\*---------------------------------------------------------------------------*/

#include "surfaceElevation.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "fvMesh.H"
#include "boundBox.H"
#include "polyBoundaryMesh.H"
#include "emptyPolyPatch.H"
#include "wedgePolyPatch.H"
#include "cyclicPolyPatch.H"
#include "processorPolyPatch.H"
#include "processorCyclicPolyPatch.H"
#include "OFstream.H"
#include "OSspecific.H"
#include "Pstream.H"
#include "PstreamReduceOps.H"
#include "DynamicList.H"
#include "reconstructionSchemes.H"
#include "IOList.H"
#include "IOdictionary.H"
#include <fstream>
#include <iomanip>

// * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(surfaceElevation, 0);
    addToRunTimeSelectionTable
    (
        functionObject,
        surfaceElevation,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * //

void Foam::functionObjects::surfaceElevation::detectBoundaryConfig()
{
    const polyBoundaryMesh& bMesh = mesh_.boundaryMesh();

    // Empty / wedge directions
    isEmpty_ = Vector<bool>(false, false, false);
    forAll(bMesh, patchi)
    {
        const polyPatch& pp = bMesh[patchi];
        if (isA<emptyPolyPatch>(pp) || isA<wedgePolyPatch>(pp))
        {
            if (pp.size() > 0)
            {
                vector nf = pp.faceAreas()[0] / mag(pp.faceAreas()[0]);
                if      (mag(nf.x()) > 0.5) isEmpty_.x() = true;
                else if (mag(nf.y()) > 0.5) isEmpty_.y() = true;
                else if (mag(nf.z()) > 0.5) isEmpty_.z() = true;
            }
        }
    }
    is2D_ = (isEmpty_.x() || isEmpty_.y() || isEmpty_.z());

    // Domain extents
    boundBox bb = mesh_.bounds();
    reduce(bb.min(), minOp<point>());
    reduce(bb.max(), maxOp<point>());

    Lx_ = bb.max().x() - bb.min().x();
    Ly_ = bb.max().y() - bb.min().y();
    Lz_ = bb.max().z() - bb.min().z();

    emptyDirThickness_ = 0;
    if      (isEmpty_.x()) emptyDirThickness_ = Lx_;
    else if (isEmpty_.y()) emptyDirThickness_ = Ly_;
    else if (isEmpty_.z()) emptyDirThickness_ = Lz_;

    // Cyclic directions (cyclic and processorCyclic counted)
    Vector<bool> cyc(false, false, false);
    forAll(bMesh, patchi)
    {
        const polyPatch& pp = bMesh[patchi];

        bool isCyc = false;
        if (isA<cyclicPolyPatch>(pp) && !isA<processorCyclicPolyPatch>(pp))
            isCyc = true;
        if (isA<processorCyclicPolyPatch>(pp))
            isCyc = true;

        if (isCyc && pp.size() > 0)
        {
            vector nf = pp.faceAreas()[0] / mag(pp.faceAreas()[0]);
            if      (mag(nf.x()) > 0.5) cyc.x() = true;
            else if (mag(nf.y()) > 0.5) cyc.y() = true;
            else if (mag(nf.z()) > 0.5) cyc.z() = true;
        }
    }
    {
        // Reduce so all procs agree (some procs may not own a cyclic patch).
        label cx = cyc.x() ? 1 : 0;
        label cy = cyc.y() ? 1 : 0;
        label cz = cyc.z() ? 1 : 0;
        reduce(cx, maxOp<label>());
        reduce(cy, maxOp<label>());
        reduce(cz, maxOp<label>());
        cyc = Vector<bool>(cx > 0, cy > 0, cz > 0);
    }
    isCyclic_ = cyc;
}


void Foam::functionObjects::surfaceElevation::markBoundaryCells(label nLayers)
{
    const polyBoundaryMesh& bMesh = mesh_.boundaryMesh();

    boundaryCells_.setSize(mesh_.nCells());
    boundaryCells_ = false;

    boolList isWallBoundaryFace(mesh_.nBoundaryFaces(), false);

    forAll(bMesh, patchi)
    {
        const polyPatch& pp = bMesh[patchi];

        if (isA<cyclicPolyPatch>(pp))      continue;
        if (isA<processorPolyPatch>(pp))   continue;
        if (isA<emptyPolyPatch>(pp))       continue;
        if (isA<wedgePolyPatch>(pp))       continue;

        forAll(pp, fi)
        {
            isWallBoundaryFace
            [
                pp.start() - mesh_.nInternalFaces() + fi
            ] = true;
        }
    }

    boolList currentLayer(mesh_.nCells(), false);
    forAll(isWallBoundaryFace, bfi)
    {
        if (isWallBoundaryFace[bfi])
        {
            const label celli =
                mesh_.faceOwner()[mesh_.nInternalFaces() + bfi];
            boundaryCells_[celli] = true;
            currentLayer[celli] = true;
        }
    }

    for (label layer = 1; layer < nLayers; layer++)
    {
        boolList nextLayer(mesh_.nCells(), false);

        for (label fi = 0; fi < mesh_.nInternalFaces(); fi++)
        {
            const label own = mesh_.faceOwner()[fi];
            const label nei = mesh_.faceNeighbour()[fi];

            if (currentLayer[own] && !boundaryCells_[nei])
            {
                boundaryCells_[nei] = true;
                nextLayer[nei] = true;
            }
            if (currentLayer[nei] && !boundaryCells_[own])
            {
                boundaryCells_[own] = true;
                nextLayer[own] = true;
            }
        }

        currentLayer = nextLayer;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::surfaceElevation::surfaceElevation
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    alphaName_("alpha.water"),
    z0_(0),
    kMax_(3),
    wetTol_(1e-6),
    dryTol_(1e-6),
    interfaceTol_(1e-12),
    statsFile_("surfaceElevation.dat"),
    writeEtaField_(true),
    gridNx_(0),
    gridNy_(0),
    gridDl_(0),
    gridXmin_(0),
    gridYmin_(0),
    gridSpecBuilt_(false),
    fineV_(0),
    isEmpty_(false, false, false),
    isCyclic_(false, false, false),
    Lx_(0), Ly_(0), Lz_(0),
    emptyDirThickness_(0),
    is2D_(false),
    boundaryCells_(mesh_.nCells(), false),
    stencilPtr_(nullptr),
    validPoints_(),
    gValid_(0),
    gIncomplete_(0),
    gCoarseLevel_(0),
    gBadBottom_(0),
    gBadTop_(0),
    gBoundarySkip_(0),
    gMultipleInterface_(0),
    etaMin_(0),
    etaMax_(0),
    etaMean_(0),
    etaVar_(0),
    etaStd_(0),
    etaSkew_(0),
    etaKurt_(0),
    etaXvar_(0),
    etaYvar_(0),
    mss_(0),
    gValidSlope_(0),
    gIncompleteSlope_(0)
{
    read(dict);

    detectBoundaryConfig();
    markBoundaryCells(kMax_);

    // amrStencil with full cubic extents (iMax = jMax = kMax = kMax_).
    // The same cubic stencil that heightFunctionAMR uses; we only read
    // the column entries topo[a3(0, 0, k)] for k = -kMax_..kMax_ and ignore
    // the lateral cells. The degenerate iMax=jMax=0 form was untested in
    // amrStencil's BFS and produced "No topology for cell" errors when the
    // column-only walk terminated before completing the stencil.
    const scalar dx = isEmpty_.x() ? emptyDirThickness_ : 0;
    const scalar dy = isEmpty_.y() ? emptyDirThickness_ : 0;
    const scalar dz = isEmpty_.z() ? emptyDirThickness_ : 0;

    stencilPtr_.reset
    (
        new amrStencil
        (
            mesh_,
            kMax_,    // iMax  -- full cubic stencil
            kMax_,    // jMax  -- full cubic stencil
            kMax_,    // kMax
            dx, dy, dz,
            isEmpty_,
            isCyclic_,
            Lx_, Ly_, Lz_
        )
    );

    Info<< "surfaceElevation: cubic stencil (2*" << kMax_ << "+1)^3"
        << "  isCyclic=(" << isCyclic_.x() << " " << isCyclic_.y()
        << " " << isCyclic_.z() << ")"
        << "  L=(" << Lx_ << " " << Ly_ << " " << Lz_ << ")" << endl;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::functionObjects::surfaceElevation::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    dict.readIfPresent("alpha",         alphaName_);
    dict.readIfPresent("z0",            z0_);
    dict.readIfPresent("kMax",          kMax_);
    dict.readIfPresent("wetTol",        wetTol_);
    dict.readIfPresent("dryTol",        dryTol_);
    dict.readIfPresent("interfaceTol",  interfaceTol_);
    dict.readIfPresent("statsFile",     statsFile_);
    dict.readIfPresent("writeEtaField", writeEtaField_);

    if (kMax_ < 1)
    {
        FatalErrorInFunction
            << "kMax must be >= 1, got " << kMax_
            << exit(FatalError);
    }

    return true;
}


bool Foam::functionObjects::surfaceElevation::execute()
{
    if (!stencilPtr_.valid())
    {
        FatalErrorInFunction
            << "amrStencil not initialised"
            << abort(FatalError);
    }

    // 1. Get alpha field
    const volScalarField& alpha =
        mesh_.lookupObject<volScalarField>(alphaName_);

    const volVectorField& C = mesh_.C();

    // 2. Use the canonical interface-cell mask from reconstructionSchemes,
    //    exactly as heightFunctionAMR does. We copy into a local boolList
    //    so the data we pass to getZoneField is byte-identical to the data
    //    we iterate over for setStencil. A const-ref directly to
    //    surf.interfaceCell() would risk invalidation if anything in the
    //    function-object machinery triggered reconstructionSchemes during
    //    the same execute() (defensive, not known to occur).
    reconstructionSchemes& surf =
        mesh_.lookupObjectRef<reconstructionSchemes>("reconstructionScheme");
    const boolList& srcMask = surf.interfaceCell();

    boolList interfaceCells(mesh_.nCells(), false);
    label nInterfaceLocal = 0;
    forAll(srcMask, celli)
    {
        if (celli < interfaceCells.size() && srcMask[celli])
        {
            interfaceCells[celli] = true;
            nInterfaceLocal++;
        }
    }

    // Diagnostic for the mask + stencil bookkeeping, useful if the run
    // aborts inside setStencil.
    {
        label gNInterface = nInterfaceLocal;
        reduce(gNInterface, sumOp<label>());
        Info<< "surfaceElevation: nInterface(global)=" << gNInterface
            << "  mesh.nCells=" << mesh_.nCells()
            << "  interfaceCells.size=" << interfaceCells.size()
            << endl;
    }

    // 3. Collect alpha values needed by the stencils on this rank
    //    (handles cross-processor traversal; populates alphaHash internally).
    Map<scalar> alphaIJK;
    const scalar t1a = mesh_.time().elapsedCpuTime();
    stencilPtr_->getZoneField(interfaceCells, alphaIJK, alpha);
    const scalar t1b = mesh_.time().elapsedCpuTime();

    // The position-keyed alpha hash. We look up alpha values by spatial
    // position rather than by cell index, mirroring heightFunctionAMR's
    // approach (the setStencil/getTopology path is not used here).
    const Map<scalar>& aHash = stencilPtr_->alphaHash();

    // Compute the finest interface-cell volume across the domain. The hash
    // is populated at the finest level only; coarse interface cells (if any
    // existed under AMR) get skipped.
    scalar fineV = GREAT;
    forAll(interfaceCells, ci)
    {
        if (interfaceCells[ci])
            fineV = min(fineV, mesh_.V()[ci]);
    }
    reduce(fineV, minOp<scalar>());
    const scalar seedVolTol = 1.5 * fineV;
    fineV_ = fineV;

    // 4. Walk interface cells and integrate the column via hash lookups.
    DynamicList<vector> validPoints;     // (x, y, etaAbs)
    validPoints.reserve(nInterfaceLocal);

    label nValid              = 0;
    label nIncomplete         = 0;
    label nCoarseLevel        = 0;
    label nBadBottom          = 0;
    label nBadTop             = 0;
    label nBoundarySkip       = 0;
    label nMultipleInterface  = 0;
    label nValidSlope         = 0;
    label nIncompleteSlope    = 0;

    // Local accumulators for the wave-statistics moments.
    scalar localSumEta    = 0;
    scalar localSumEta2   = 0;
    scalar localSumEta3   = 0;
    scalar localSumEta4   = 0;
    scalar localEtaMin    =  GREAT;
    scalar localEtaMax    = -GREAT;

    // Slope variance accumulators
    scalar localSumEtaX2  = 0;
    scalar localSumEtaY2  = 0;

    const scalar t2a = mesh_.time().elapsedCpuTime();

    // Helper: compute column eta at a horizontal offset (dx,dy) from seedCC
    // by walking the 7-cell column via alphaHash lookups. Returns true if all
    // 7 lookups succeeded. The integrity check (wet bottom, dry top) is NOT
    // applied to lateral columns; their purpose is the central difference.
    auto columnEta = [&]
    (
        const point& centerCC,
        scalar dl,
        scalar dx,
        scalar dy,
        scalar& etaOut
    ) -> bool
    {
        scalar sum = 0;
        for (label k = -kMax_; k <= kMax_; k++)
        {
            point eCC = centerCC;
            eCC.x() += dx;
            eCC.y() += dy;
            eCC.z() += k * dl;
            eCC = stencilPtr_->cyclicWrap(eCC);

            const label key = stencilPtr_->posKey(eCC);
            const auto iter = aHash.cfind(key);
            if (!iter.good()) return false;
            sum += iter.val();
        }

        const scalar zBotFace = centerCC.z() - (scalar(kMax_) + 0.5) * dl;
        etaOut = zBotFace + sum * dl;
        return true;
    };

    forAll(interfaceCells, celli)
    {
        if (!interfaceCells[celli]) continue;

        // Skip cells too close to a non-cyclic wall: column would leave domain
        if (celli < boundaryCells_.size() && boundaryCells_[celli])
        {
            nBoundarySkip++;
            continue;
        }

        // Skip non-finest cells (alphaHash is built at finest level only)
        if (mesh_.V()[celli] > seedVolTol)
        {
            nCoarseLevel++;
            continue;
        }

        const scalar V = mesh_.V()[celli];
        scalar dl;
        if (is2D_)
            dl = Foam::sqrt(V / emptyDirThickness_);
        else
            dl = Foam::cbrt(V);

        const point seedCC = mesh_.C()[celli];

        // Read all 7 column alpha values via spatial-hash lookup.
        List<scalar> A(2*kMax_ + 1, 0.0);
        bool incomplete = false;

        for (label k = -kMax_; k <= kMax_ && !incomplete; k++)
        {
            point eCC = seedCC;
            eCC.z() += k * dl;
            eCC = stencilPtr_->cyclicWrap(eCC);

            const label key = stencilPtr_->posKey(eCC);
            const auto iter = aHash.cfind(key);

            if (iter.good())
                A[k + kMax_] = iter.val();
            else
                incomplete = true;
        }

        if (incomplete)
        {
            nIncomplete++;
            continue;
        }

        // Integrity: bottom must be fully wet, top must be fully dry
        const scalar aBot = A[0];
        const scalar aTop = A[2*kMax_];

        if (aBot < 1.0 - wetTol_)
        {
            nBadBottom++;
            if (aBot > interfaceTol_ && aBot < 1.0 - interfaceTol_)
                nMultipleInterface++;
            continue;
        }
        if (aTop > dryTol_)
        {
            nBadTop++;
            if (aTop > interfaceTol_ && aTop < 1.0 - interfaceTol_)
                nMultipleInterface++;
            continue;
        }

        // Column sum at the seed
        scalar etaLocal = 0;
        for (label k = 0; k < 2*kMax_ + 1; k++)
            etaLocal += A[k];
        etaLocal *= dl;

        const scalar zBotFace = seedCC.z() - (scalar(kMax_) + 0.5) * dl;
        const scalar etaAbs       = zBotFace + etaLocal;
        const scalar etaExcursion = etaAbs - z0_;

        validPoints.append(vector(seedCC.x(), seedCC.y(), etaAbs));

        const scalar e2 = etaExcursion * etaExcursion;
        localSumEta  += etaExcursion;
        localSumEta2 += e2;
        localSumEta3 += e2 * etaExcursion;
        localSumEta4 += e2 * e2;
        if (etaExcursion < localEtaMin) localEtaMin = etaExcursion;
        if (etaExcursion > localEtaMax) localEtaMax = etaExcursion;

        nValid++;

        // ----- Slope statistics: lateral columns at (+/-dl, 0) and (0, +/-dl)
        // If any lateral column has a missing hash entry, drop the seed from
        // slope statistics (but keep it in eta statistics above). Slopes
        // require the four lateral neighbours; etaXY skipped if non-3D.
        if (is2D_)
        {
            // 2D: only one lateral direction is meaningful; skip slopes
            continue;
        }

        scalar etaXp, etaXm, etaYp, etaYm;
        const bool slopeOK =
            columnEta(seedCC, dl, +dl,   0, etaXp) &&
            columnEta(seedCC, dl, -dl,   0, etaXm) &&
            columnEta(seedCC, dl,   0, +dl, etaYp) &&
            columnEta(seedCC, dl,   0, -dl, etaYm);

        if (!slopeOK)
        {
            nIncompleteSlope++;
            continue;
        }

        const scalar etaX = (etaXp - etaXm) / (2.0 * dl);
        const scalar etaY = (etaYp - etaYm) / (2.0 * dl);
        localSumEtaX2 += etaX * etaX;
        localSumEtaY2 += etaY * etaY;
        nValidSlope++;
    }

    const scalar t2b = mesh_.time().elapsedCpuTime();


    // 5. Reduce diagnostics
    label gValid             = nValid;
    label gIncomplete        = nIncomplete;
    label gCoarseLevel       = nCoarseLevel;
    label gBadBottom         = nBadBottom;
    label gBadTop            = nBadTop;
    label gBoundarySkip      = nBoundarySkip;
    label gMultipleInterface = nMultipleInterface;
    label gValidSlope        = nValidSlope;
    label gIncompleteSlope   = nIncompleteSlope;

    reduce(gValid,             sumOp<label>());
    reduce(gIncomplete,        sumOp<label>());
    reduce(gCoarseLevel,       sumOp<label>());
    reduce(gBadBottom,         sumOp<label>());
    reduce(gBadTop,            sumOp<label>());
    reduce(gBoundarySkip,      sumOp<label>());
    reduce(gMultipleInterface, sumOp<label>());
    reduce(gValidSlope,        sumOp<label>());
    reduce(gIncompleteSlope,   sumOp<label>());

    scalar gSumEta  = localSumEta;
    scalar gSumEta2 = localSumEta2;
    scalar gSumEta3 = localSumEta3;
    scalar gSumEta4 = localSumEta4;
    scalar gEtaMin  = localEtaMin;
    scalar gEtaMax  = localEtaMax;
    reduce(gSumEta,  sumOp<scalar>());
    reduce(gSumEta2, sumOp<scalar>());
    reduce(gSumEta3, sumOp<scalar>());
    reduce(gSumEta4, sumOp<scalar>());
    reduce(gEtaMin,  minOp<scalar>());
    reduce(gEtaMax,  maxOp<scalar>());

    // Slope-sum reductions
    scalar gSumEtaX2 = localSumEtaX2;
    scalar gSumEtaY2 = localSumEtaY2;
    reduce(gSumEtaX2, sumOp<scalar>());
    reduce(gSumEtaY2, sumOp<scalar>());

    // Statistical moments. Normalised by nValid (uniform-cell average).
    //   m1 = <eta>     = sum_eta / N
    //   m2 = <eta^2>   = sum_eta2 / N      etc.
    // Central moments via the binomial expansion:
    //   var  = m2 - m1^2
    //   mu3  = m3 - 3*m1*m2 + 2*m1^3
    //   mu4  = m4 - 4*m1*m3 + 6*m1^2*m2 - 3*m1^4
    //   skew = mu3 / std^3
    //   kurt = mu4 / std^4    (NOT excess; subtract 3 in Python if you prefer)
    scalar etaMean = 0, etaVar = 0, etaStd = 0, etaSkew = 0, etaKurt = 0;
    if (gValid > 0)
    {
        const scalar N  = scalar(gValid);
        const scalar m1 = gSumEta  / N;
        const scalar m2 = gSumEta2 / N;
        const scalar m3 = gSumEta3 / N;
        const scalar m4 = gSumEta4 / N;
        etaMean = m1;
        etaVar  = max(0.0, m2 - m1*m1);
        etaStd  = Foam::sqrt(etaVar);
        if (etaStd > SMALL)
        {
            const scalar mu3 = m3 - 3*m1*m2 + 2*m1*m1*m1;
            const scalar mu4 = m4 - 4*m1*m3 + 6*m1*m1*m2 - 3*m1*m1*m1*m1;
            etaSkew = mu3 / (etaStd * etaStd * etaStd);
            etaKurt = mu4 / (etaVar * etaVar);
        }
    }

    // Slope variances. <eta_x^2> and <eta_y^2> are accumulated unweighted
    // over cells where all four lateral columns succeeded. mss is the total
    // mean square slope <|grad eta|^2>.
    scalar etaXvar = 0, etaYvar = 0, mss = 0;
    if (gValidSlope > 0)
    {
        const scalar Ns = scalar(gValidSlope);
        etaXvar = gSumEtaX2 / Ns;
        etaYvar = gSumEtaY2 / Ns;
        mss     = etaXvar + etaYvar;
    }

    // 6. Stash results in member state for write() to flush
    validPoints_.transfer(validPoints);
    gValid_             = gValid;
    gIncomplete_        = gIncomplete;
    gCoarseLevel_       = gCoarseLevel;
    gBadBottom_         = gBadBottom;
    gBadTop_            = gBadTop;
    gBoundarySkip_      = gBoundarySkip;
    gMultipleInterface_ = gMultipleInterface;
    etaMin_             = gEtaMin;
    etaMax_             = gEtaMax;
    etaMean_            = etaMean;
    etaVar_             = etaVar;
    etaStd_             = etaStd;
    etaSkew_            = etaSkew;
    etaKurt_            = etaKurt;
    etaXvar_            = etaXvar;
    etaYvar_            = etaYvar;
    mss_                = mss;
    gValidSlope_        = gValidSlope;
    gIncompleteSlope_   = gIncompleteSlope;

    // 7. Console diagnostics
    scalar dt1 = t1b - t1a;
    scalar dt2 = t2b - t2a;
    reduce(dt1, maxOp<scalar>());
    reduce(dt2, maxOp<scalar>());

    Info<< "surfaceElevation: t=" << mesh_.time().timeName()
        << "  valid=" << gValid
        << "  fail=("
        <<   "incomp="     << gIncomplete
        << " coarseLvl="   << gCoarseLevel
        << " badBot="      << gBadBottom
        << " badTop="      << gBadTop
        << " bndry="       << gBoundarySkip
        << " multi="       << gMultipleInterface
        << ")"
        << "  etaMin="     << gEtaMin
        << "  etaMax="     << gEtaMax
        << "  etaStd="     << etaStd
        << "  mss="        << mss
        << "  T1(getZ)="   << dt1
        << "s T2(loop)="   << dt2 << "s"
        << endl;

    // Build the binned eta grid on master and broadcast to all ranks.
    // Done every execute so downstream consumers (zetaSampler) see a
    // current eta on their own cadence.
    buildEtaGrid();

    return true;
}


// Write the column header to the time-series .dat file.
void Foam::functionObjects::surfaceElevation::writeStatsHeader
(
    std::ofstream& os
) const
{
    os << "# time"
       << " sigmaEta Hs etaMin etaMax HmaxOverHs skew kurt"
       << " etaXvar etaYvar mss"
       << " nValid nValidSlope nIncomp nIncompSlope nCoarseLvl"
       << " nBadBot nBadTop nBndry nMulti"
       << "\n";
}


// Write the regular-grid spec once at startup. Builds the grid bounds
// from validPoints_ on master. The spec lives at
//   <case>/postProcessing/surfaceElevation/grid.txt
void Foam::functionObjects::surfaceElevation::writeGridSpec() const
{
    if (!Pstream::master()) return;
    if (gridSpecBuilt_) return;

    const fileName gridDir
    (
        mesh_.time().globalPath() / "postProcessing" / name()
    );
    mkDir(gridDir);

    std::ofstream os((gridDir / "grid.txt").c_str());
    os << "# surfaceElevation grid spec\n"
       << "# eta.bin layout: " << gridNx_ << " x " << gridNy_
       << " doubles, row-major, col-major in (i,j) order\n"
       << "# i = round((x - xmin) / dl)\n"
       << "# j = round((y - ymin) / dl)\n"
       << "# Missing cells are written as NaN.\n"
       << "Nx   " << gridNx_   << "\n"
       << "Ny   " << gridNy_   << "\n"
       << "dl   " << std::scientific << std::setprecision(12)
                  << gridDl_   << "\n"
       << "xmin " << gridXmin_ << "\n"
       << "ymin " << gridYmin_ << "\n";
    os.close();

    gridSpecBuilt_ = true;
}


// Gather validPoints_ from all ranks to master, bin onto the regular Lvl1
// grid, write to <case>/postProcessing/surfaceElevation/<time>/eta.bin.
// Master-only single binary file; missing cells are NaN.
void Foam::functionObjects::surfaceElevation::buildEtaGrid()
{
    // Called every execute(). Master gathers validPoints_ from all ranks,
    // bins onto regular Nx x Ny grid, broadcasts the resulting grid back
    // to all ranks, and stores it in a registered IOList<scalar> on the
    // mesh objectRegistry. Downstream functionObjects pick it up by name
    // via mesh.lookupObject<IOList<scalar>>(etaGridName()). The grid
    // spec (Nx, Ny, dl, xmin, ymin) lives in a parallel IOdictionary
    // under etaGridSpecName(). Both are NO_WRITE.
    //
    // Also updates the volScalarField "zeta" = z - eta(x, y, t) on every
    // cell. The zeta field is registered to the mesh once, NO_WRITE; the
    // user can sample it via standard surfaces / isoSurface fOs.

    // Gather validPoints_ from all ranks to master
    List<List<vector>> allPoints(Pstream::nProcs());
    allPoints[Pstream::myProcNo()] = validPoints_;
    Pstream::gatherList(allPoints);

    if (Pstream::master())
    {
        // Build grid spec on first call
        if (!gridSpecBuilt_)
        {
            scalar xmin = GREAT, xmax = -GREAT;
            scalar ymin = GREAT, ymax = -GREAT;

            forAll(allPoints, p)
            {
                const List<vector>& pts = allPoints[p];
                forAll(pts, i)
                {
                    const scalar x = pts[i].x();
                    const scalar y = pts[i].y();
                    xmin = min(xmin, x);  xmax = max(xmax, x);
                    ymin = min(ymin, y);  ymax = max(ymax, y);
                }
            }

            scalar dl;
            if (is2D_)
                dl = Foam::sqrt(fineV_ / emptyDirThickness_);
            else
                dl = Foam::cbrt(fineV_);

            gridDl_   = dl;
            gridXmin_ = xmin;
            gridYmin_ = ymin;
            gridNx_   = label(round((xmax - xmin) / dl)) + 1;
            gridNy_   = label(round((ymax - ymin) / dl)) + 1;

            Info<< "surfaceElevation: grid spec  Nx=" << gridNx_
                << " Ny=" << gridNy_
                << " dl="   << gridDl_
                << " xmin=" << gridXmin_
                << " ymin=" << gridYmin_ << endl;

            writeGridSpec();
        }
    }

    // Broadcast spec to all ranks
    Pstream::broadcast(gridNx_);
    Pstream::broadcast(gridNy_);
    Pstream::broadcast(gridDl_);
    Pstream::broadcast(gridXmin_);
    Pstream::broadcast(gridYmin_);

    // Lazy-construct the registered IOList for the eta grid, the
    // IOdictionary for the spec, and the volScalarField zeta. All are
    // NO_WRITE so they never go to disk through the standard mechanism.
    if (!etaGridIOPtr_)
    {
        etaGridIOPtr_.reset
        (
            new IOList<scalar>
            (
                IOobject
                (
                    etaGridName(),
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    IOobject::REGISTER
                ),
                List<scalar>(gridNx_ * gridNy_,
                             std::numeric_limits<scalar>::quiet_NaN())
            )
        );
    }
    if (!etaGridSpecPtr_)
    {
        etaGridSpecPtr_.reset
        (
            new IOdictionary
            (
                IOobject
                (
                    etaGridSpecName(),
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    IOobject::REGISTER
                )
            )
        );
        etaGridSpecPtr_->add("Nx",   gridNx_);
        etaGridSpecPtr_->add("Ny",   gridNy_);
        etaGridSpecPtr_->add("dl",   gridDl_);
        etaGridSpecPtr_->add("xmin", gridXmin_);
        etaGridSpecPtr_->add("ymin", gridYmin_);
    }
    if (!zetaPtr_)
    {
        zetaPtr_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    zetaName(),
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    IOobject::REGISTER
                ),
                mesh_,
                dimensionedScalar("zeta", dimLength, scalar(0))
            )
        );
    }

    // Resize and bin into the registered IOList in place
    etaGridIOPtr_->resize_nocopy(gridNx_ * gridNy_);
    List<scalar>& etaArr = *etaGridIOPtr_;
    if (Pstream::master())
    {
        const scalar nanVal = std::numeric_limits<scalar>::quiet_NaN();
        for (label i = 0; i < etaArr.size(); ++i) etaArr[i] = nanVal;

        forAll(allPoints, p)
        {
            const List<vector>& pts = allPoints[p];
            forAll(pts, i)
            {
                const label ix =
                    label(round((pts[i].x() - gridXmin_) / gridDl_));
                const label iy =
                    label(round((pts[i].y() - gridYmin_) / gridDl_));
                if (ix < 0 || ix >= gridNx_ || iy < 0 || iy >= gridNy_)
                    continue;
                etaArr[iy * gridNx_ + ix] = pts[i].z() - z0_;
            }
        }
    }

    // Broadcast the grid array to all ranks
    Pstream::broadcast(etaArr);

    // Populate zeta = z - eta(x_c, y_c) at every cell using the etaArr
    // lookup. Cells outside the binned (x, y) extent get zero (rare;
    // grid was built from the extent of validPoints_).
    {
        volScalarField& zeta = *zetaPtr_;
        const vectorField& C = mesh_.C().internalField();
        scalarField& zetaI = zeta.primitiveFieldRef();

        for (label celli = 0; celli < zetaI.size(); ++celli)
        {
            const scalar xc = C[celli].x();
            const scalar yc = C[celli].y();
            const label ix = label(round((xc - gridXmin_) / gridDl_));
            const label iy = label(round((yc - gridYmin_) / gridDl_));

            scalar etaLocal = 0;
            if (ix >= 0 && ix < gridNx_ && iy >= 0 && iy < gridNy_)
            {
                const scalar v = etaArr[iy * gridNx_ + ix];
                if (!std::isnan(v)) etaLocal = v;
            }
            zetaI[celli] = C[celli].z() - (z0_ + etaLocal);
        }
        zeta.correctBoundaryConditions();
    }

    gridSpecBuilt_ = true;
}


void Foam::functionObjects::surfaceElevation::writeEtaBinary() const
{
    if (!Pstream::master()) return;
    if (!gridSpecBuilt_) return;
    if (!etaGridIOPtr_) return;

    const fileName outDir
    (
        mesh_.time().globalPath() / "postProcessing" /
        name() / mesh_.time().timeName()
    );
    mkDir(outDir);

    const List<scalar>& etaArr = *etaGridIOPtr_;
    std::ofstream os((outDir / "eta.bin").c_str(), std::ios::binary);
    os.write
    (
        reinterpret_cast<const char*>(etaArr.cdata()),
        std::streamsize(etaArr.size() * sizeof(scalar))
    );
    os.close();

    Info<< "surfaceElevation: wrote eta.bin"
        << "  size=" << (etaArr.size()*sizeof(scalar))/1024/1024 << "MB"
        << endl;
}


bool Foam::functionObjects::surfaceElevation::write()
{
    // Skip if execute() has not produced any data yet.
    if (!stencilPtr_.valid())
    {
        return true;
    }

    const scalar t3a = mesh_.time().elapsedCpuTime();

    // 1. Master-only: append a row to the wave-statistics file at case root.
    //    Always append, never truncate. Restart picks up where it left
    //    off; the user can manually rotate the file between runs if
    //    they want clean separation.
    if (Pstream::master())
    {
        const fileName statsPath
        (
            mesh_.time().globalPath() / statsFile_
        );

        // Header is written only if the file is empty (or absent).
        std::ifstream check(statsPath.c_str(), std::ios::ate);
        const bool isEmpty = !check.good() || check.tellg() == 0;
        check.close();

        std::ofstream os(statsPath.c_str(), std::ios::out | std::ios::app);
        if (isEmpty)
        {
            writeStatsHeader(os);
        }

        const scalar Hs         = 4.0 * etaStd_;
        const scalar Hmax       = etaMax_ - etaMin_;
        const scalar HmaxOverHs = (Hs > SMALL) ? Hmax / Hs : 0.0;

        os << std::scientific << std::setprecision(8)
           << mesh_.time().value() << " "
           << etaStd_   << " "
           << Hs        << " "
           << etaMin_   << " "
           << etaMax_   << " "
           << HmaxOverHs << " "
           << etaSkew_  << " "
           << etaKurt_  << " "
           << etaXvar_  << " "
           << etaYvar_  << " "
           << mss_      << " "
           << gValid_             << " "
           << gValidSlope_        << " "
           << gIncomplete_        << " "
           << gIncompleteSlope_   << " "
           << gCoarseLevel_       << " "
           << gBadBottom_         << " "
           << gBadTop_            << " "
           << gBoundarySkip_      << " "
           << gMultipleInterface_ << "\n";
        os.close();
    }

    // 2. Optional: dump binary eta(x, y) snapshot. Master-only single file.
    if (writeEtaField_)
    {
        writeEtaBinary();
    }

    const scalar t3b = mesh_.time().elapsedCpuTime();
    scalar dt3 = t3b - t3a;
    reduce(dt3, maxOp<scalar>());

    Info<< "surfaceElevation: wrote t=" << mesh_.time().timeName()
        << "  T3(IO)=" << dt3 << "s" << endl;

    return true;
}


// ************************************************************************* //
