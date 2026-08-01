/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2026 Asim Onder
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of the unofficial extension library geometricVoFCartesian.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "mycCartesianAMR.H"
#include "fvc.H"
#include "addToRunTimeSelectionTable.H"
#include "emptyPolyPatch.H"
#include "wedgePolyPatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace reconstruction
{
    defineTypeNameAndDebug(mycCartesianAMR, 0);
    addToRunTimeSelectionTable
    (
        reconstructionSchemes,
        mycCartesianAMR,
        components
    );
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::reconstruction::mycCartesianAMR::markBoundaryCells(label nLayers)
{
    const polyBoundaryMesh& bMesh = mesh_.boundaryMesh();
    boolList isWallBoundaryFace(mesh_.nBoundaryFaces(), false);
    forAll(bMesh, patchi)
    {
        const polyPatch& pp = bMesh[patchi];
        if (isA<cyclicPolyPatch>(pp)) continue;
        if (isA<processorPolyPatch>(pp)) continue;
        if (isA<emptyPolyPatch>(pp)) continue;
        if (isA<wedgePolyPatch>(pp)) continue;
        forAll(pp, fi)
            isWallBoundaryFace[pp.start()-mesh_.nInternalFaces()+fi] = true;
    }
    boolList currentLayer(mesh_.nCells(), false);
    forAll(isWallBoundaryFace, bfi)
    {
        if (isWallBoundaryFace[bfi])
        {
            const label celli = mesh_.faceOwner()[mesh_.nInternalFaces()+bfi];
            boundaryCells_[celli] = true;
            currentLayer[celli] = true;
        }
    }
    for (label layer=1; layer<nLayers; layer++)
    {
        boolList nextLayer(mesh_.nCells(), false);
        for (label fi=0; fi<mesh_.nInternalFaces(); fi++)
        {
            const label o=mesh_.faceOwner()[fi], n=mesh_.faceNeighbour()[fi];
            if (currentLayer[o] && !boundaryCells_[n])
            { boundaryCells_[n]=true; nextLayer[n]=true; }
            if (currentLayer[n] && !boundaryCells_[o])
            { boundaryCells_[o]=true; nextLayer[o]=true; }
        }
        currentLayer = nextLayer;
    }
}


void Foam::reconstruction::mycCartesianAMR::gradSurf(const volScalarField& phi)
{
    Map<scalar> phiIJK;
    amrStencilPtr_->getZoneField(interfaceCell_, phiIJK, phi);

    // gradient normals as the fallback for cells whose cube the
    // spatial hash cannot complete (band edges, unresolved remotes)
    const volVectorField gradAlpha(fvc::grad(phi));

    if (debugCheck_)
    {
        // findCell lazily builds tetBasePtIs, whose construction
        // synchronizes across ranks; build it collectively here so
        // per-cell calls cannot deadlock ranks that own no
        // interface cells
        (void)mesh_.tetBasePtIs();
    }

    const bool twoD = isEmpty_.x() || isEmpty_.y() || isEmpty_.z();

    label youngCount = 0;
    label CCDCount = 0;
    label fallbackCount = 0;

    // diagnostic: cube slots not resolvable by local findCell, i.e.
    // relay-filled entries, queued here and verified against their
    // owning rank after the loop
    DynamicList<point> dbgQPts;
    DynamicList<scalar> dbgQVals;

    forAll(interfaceLabels_, i)
    {
        const label celli = interfaceLabels_[i];

        if (boundaryCells_[celli])
        {
            scalar alphaNew = (alpha1_[celli] < 0.5) ? 0.0 : 1.0;
            Pout<< "mycCartesianAMR: cell at " << mesh_.C()[celli]
                << " with alpha=" << alpha1_[celli]
                << " is at a non-cyclic boundary; setting alpha to "
                << alphaNew << nl;
            alpha1_[celli] = alphaNew;
            interfaceNormal_[i] = vector::zero;
            continue;
        }

        // local fine spacing; in 2D the empty-direction thickness must
        // be divided out of the cell volume
        const scalar h =
            twoD
          ? Foam::sqrt(mesh_.V()[celli]/emptyDirThickness_)
          : Foam::cbrt(mesh_.V()[celli]);

        const label nMiss =
            amrStencilPtr_->setStencilFromHash(mesh_.C()[celli], h);

        // temporary diagnostic (amrDebugCheck in the alpha dict):
        // verify every locally resolvable cube slot against the
        // mesh's own findCell answer, catching wrong-value hash
        // entries and misaligned bins
        if (debugCheck_)
        {
            static label dbgBad = 0;
            const List<scalar> Acube(amrStencilPtr_->getStencil());
            for (label kk = -1; kk <= 1; ++kk)
            {
                if (isEmpty_.z() && kk != 0) continue;
                for (label jj = -1; jj <= 1; ++jj)
                {
                    if (isEmpty_.y() && jj != 0) continue;
                    for (label ii = -1; ii <= 1; ++ii)
                    {
                        if (isEmpty_.x() && ii != 0) continue;
                        point p = mesh_.C()[celli]
                            + vector(ii*h, jj*h, kk*h);
                        p = amrStencilPtr_->cyclicWrap(p);
                        const label c2 = mesh_.findCell(p);
                        const label idx = amrStencilPtr_->a3(ii, jj, kk);
                        if (c2 < 0)
                        {
                            // relay-filled slot: verify remotely
                            dbgQPts.append(p);
                            dbgQVals.append(Acube[idx]);
                            continue;
                        }
                        if
                        (
                            mag(Acube[idx] - alpha1_[c2]) > 1e-10
                         && dbgBad < 40
                        )
                        {
                            ++dbgBad;
                            Pout<< "mycCartesianAMR BADSLOT C="
                                << mesh_.C()[celli]
                                << " off=(" << ii << " " << jj
                                << " " << kk << ")"
                                << " hash=" << Acube[idx]
                                << " true=" << alpha1_[c2] << nl;
                        }
                    }
                }
            }
        }

        if (nMiss > 0)
        {
            // temporary diagnostic: report where the incomplete cubes
            // sit, to attribute them to proc corners, the cyclic seam,
            // or hash misalignment
            static label dbgReports = 0;
            if (dbgReports < 30)
            {
                ++dbgReports;
                Pout<< "mycCartesianAMR miss: C=" << mesh_.C()[celli]
                    << " nMiss=" << nMiss << " h=" << h
                    << " a=" << alpha1_[celli] << nl;
            }
            vector g = gradAlpha[celli];
            const scalar mg1 = fabs(g.x()) + fabs(g.y()) + fabs(g.z());
            interfaceNormal_[i] = (mg1 > 0) ? vector(g/mg1) : vector::zero;
            fallbackCount += 1;
            continue;
        }

        vector m1 = amrStencilPtr_->calcYoungNormal();
        vector m2 = amrStencilPtr_->calcCCDNormal();
        const vector m2C = m2;

        if (Foam::mag(m1) + Foam::mag(m2) == 0)
        {
            scalar alphaNew = (alpha1_[celli] < 0.5) ? 0.0 : 1.0;
            Pout<< "mycCartesianAMR: Young and CCD normals both zero at "
                << mesh_.C()[celli] << " with alpha=" << alpha1_[celli]
                << "; setting alpha to " << alphaNew << nl;
            alpha1_[celli] = alphaNew;
            interfaceNormal_[i] = vector::zero;
            continue;
        }

        if (Foam::mag(m1) > 0)
        {
            m1 /= (fabs(m1.x()) + fabs(m1.y()) + fabs(m1.z()));
        }
        if (Foam::mag(m2) > 0)
        {
            m2 /= (fabs(m2.x()) + fabs(m2.y()) + fabs(m2.z()));
        }

        if (Foam::mag(m2C) == 0)
        {
            interfaceNormal_[i] = m1;
            youngCount += 1;
        }
        else if (fabs(m2C.y()) == 1.0)
        {
            if (fabs(m1.y()) < fabs(m2.y()))
            {
                interfaceNormal_[i] = m1;
                youngCount += 1;
            }
            else
            {
                interfaceNormal_[i] = m2;
                CCDCount += 1;
            }
        }
        else if (fabs(m2C.z()) == 1.0)
        {
            if (fabs(m1.z()) < fabs(m2.z()))
            {
                interfaceNormal_[i] = m1;
                youngCount += 1;
            }
            else
            {
                interfaceNormal_[i] = m2;
                CCDCount += 1;
            }
        }
        else if (fabs(m2C.x()) == 1.0)
        {
            if (fabs(m1.x()) < fabs(m2.x()))
            {
                interfaceNormal_[i] = m1;
                youngCount += 1;
            }
            else
            {
                interfaceNormal_[i] = m2;
                CCDCount += 1;
            }
        }
        else
        {
            // incomplete cube or degenerate CCD state under AMR:
            // degrade to the Youngs normal instead of aborting
            // (mycCartesian aborts here); zero-normal snap only when
            // Youngs is also unavailable
            if (Foam::mag(m1) > 0)
            {
                interfaceNormal_[i] = m1;
                fallbackCount += 1;
            }
            else
            {
                interfaceNormal_[i] = m2;
                fallbackCount += 1;
            }
        }
    }

    // diagnostic: every relay-filled slot goes back to the rank that
    // owns the queried position and is compared against its field;
    // buffered all-to-all, the library's own proven exchange pattern
    if (debugCheck_ && Pstream::parRun())
    {
        const pointField myPts(dbgQPts);
        const scalarField myVals(dbgQVals);

        PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);
        for (const int proci : Pstream::allProcs())
        {
            if (proci == Pstream::myProcNo()) continue;
            UOPstream toProc(proci, pBufs);
            toProc << myPts << myVals;
        }
        pBufs.finishedSends();

        static label dbgRemBad = 0;
        label nChecked = 0;
        for (const int proci : Pstream::allProcs())
        {
            if (proci == Pstream::myProcNo()) continue;
            UIPstream fromProc(proci, pBufs);
            pointField pts(fromProc);
            scalarField vals(fromProc);
            forAll(pts, qi)
            {
                const label c2 = mesh_.findCell(pts[qi]);
                if (c2 < 0) continue;
                ++nChecked;
                if
                (
                    mag(vals[qi] - alpha1_[c2]) > 1e-10
                 && dbgRemBad < 40
                )
                {
                    ++dbgRemBad;
                    Pout<< "mycCartesianAMR BADREMOTE p=" << pts[qi]
                        << " hash=" << vals[qi]
                        << " true=" << alpha1_[c2] << nl;
                }
            }
        }
        reduce(nChecked, sumOp<label>());
        Info<< "mycCartesianAMR: remote slots checked = "
            << nChecked << endl;
    }

    if (Pstream::parRun())
    {
        Foam::reduce(youngCount, sumOp<label>());
        Foam::reduce(CCDCount, sumOp<label>());
        Foam::reduce(fallbackCount, sumOp<label>());
    }

    Info<< "mycCartesianAMR: youngCount=" << youngCount
        << " CCDCount=" << CCDCount
        << " fallbackCount=" << fallbackCount << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::reconstruction::mycCartesianAMR::mycCartesianAMR
(
    volScalarField& alpha1,
    const surfaceScalarField& phi,
    const volVectorField& U,
    const dictionary& dict
)
:
    reconstructionSchemes
    (
        typeName,
        alpha1,
        phi,
        U,
        dict
    ),
    mesh_(alpha1.mesh()),
    interfaceNormal_(fvc::grad(alpha1)),
    isoFaceTol_(modelDict().lookupOrDefault<scalar>("isoFaceTol", 1e-8)),
    surfCellTol_(modelDict().lookupOrDefault<scalar>("surfCellTol", 1e-8)),
    sIterPLIC_(mesh_, surfCellTol_),
    isEmpty_(false, false, false),
    emptyDirThickness_(0),
    boundaryCells_(mesh_.nCells(), false),
    debugCheck_(modelDict().lookupOrDefault<bool>("amrDebugCheck", false))
{
    Info<< "Constructing the interface using mycCartesianAMR method..."
        << endl;

    // empty (2D) directions
    const polyBoundaryMesh& bMesh = mesh_.boundaryMesh();
    forAll(bMesh, patchi)
    {
        const polyPatch& pp = bMesh[patchi];
        if (isA<emptyPolyPatch>(pp) || isA<wedgePolyPatch>(pp))
        {
            if (pp.size() > 0)
            {
                vector nf = pp.faceAreas()[0]/mag(pp.faceAreas()[0]);
                if (mag(nf.x()) > 0.5) isEmpty_.x() = true;
                else if (mag(nf.y()) > 0.5) isEmpty_.y() = true;
                else if (mag(nf.z()) > 0.5) isEmpty_.z() = true;
            }
        }
    }

    boundBox bb = mesh_.bounds();
    reduce(bb.min(), minOp<point>());
    reduce(bb.max(), maxOp<point>());
    if (isEmpty_.x()) emptyDirThickness_ = bb.max().x() - bb.min().x();
    else if (isEmpty_.y()) emptyDirThickness_ = bb.max().y() - bb.min().y();
    else if (isEmpty_.z()) emptyDirThickness_ = bb.max().z() - bb.min().z();

    // cyclic directions
    Vector<bool> isCyclic(false, false, false);
    forAll(bMesh, patchi)
    {
        const polyPatch& pp = bMesh[patchi];
        bool isCyclicPatch = false;
        if (isA<cyclicPolyPatch>(pp) && !isA<processorCyclicPolyPatch>(pp))
            isCyclicPatch = true;
        if (isA<processorCyclicPolyPatch>(pp))
            isCyclicPatch = true;
        if (isCyclicPatch)
        {
            const vectorField& fAreas = pp.faceAreas();
            if (fAreas.size() > 0)
            {
                vector nf = fAreas[0]/mag(fAreas[0]);
                if (mag(nf.x()) > 0.5) isCyclic.x() = true;
                else if (mag(nf.y()) > 0.5) isCyclic.y() = true;
                else if (mag(nf.z()) > 0.5) isCyclic.z() = true;
            }
        }
    }
    {
        label cx=isCyclic.x()?1:0, cy=isCyclic.y()?1:0, cz=isCyclic.z()?1:0;
        reduce(cx, maxOp<label>());
        reduce(cy, maxOp<label>());
        reduce(cz, maxOp<label>());
        isCyclic.x()=(cx>0); isCyclic.y()=(cy>0); isCyclic.z()=(cz>0);
    }

    const scalar Lx = bb.max().x() - bb.min().x();
    const scalar Ly = bb.max().y() - bb.min().y();
    const scalar Lz = bb.max().z() - bb.min().z();

    Info<< "mycCartesianAMR: isCyclic=("
        << isCyclic.x() << " " << isCyclic.y() << " " << isCyclic.z()
        << ") L=(" << Lx << " " << Ly << " " << Lz << ")" << endl;

    // The first spacing argument becomes uniformStencil::dx_. The
    // CCD normal mixes literal unity components with slope terms
    // divided by dx_, so dx_ must be the TRUE fine spacing of the
    // interface band; an arbitrary scale skews the Youngs-CCD
    // selection and the normals themselves. The band sits at the
    // finest refinement level, whose spacing is constant through a
    // run, so the minimum cell size at construction is exact.
    scalar hFine = GREAT;
    {
        const bool twoD0 = isEmpty_.x() || isEmpty_.y() || isEmpty_.z();
        forAll(mesh_.V(), celli)
        {
            const scalar cs =
                twoD0
              ? Foam::sqrt(mesh_.V()[celli]/emptyDirThickness_)
              : Foam::cbrt(mesh_.V()[celli]);
            hFine = min(hFine, cs);
        }
        reduce(hFine, minOp<scalar>());
    }
    Info<< "mycCartesianAMR: fine spacing dx = " << hFine << endl;

    const scalar dy = isEmpty_.y() ? emptyDirThickness_ : 0;
    const scalar dz = isEmpty_.z() ? emptyDirThickness_ : 0;
    amrStencilPtr_.reset
    (
        new amrStencil(mesh_, 1, 1, 1,
            hFine, dy, dz, isEmpty_, isCyclic, Lx, Ly, Lz)
    );

    markBoundaryCells(1);

    reconstruct();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::reconstruction::mycCartesianAMR::reconstruct(bool forceUpdate)
{
    const bool uptodate = alreadyReconstructed(forceUpdate);

    if (uptodate && !forceUpdate)
    {
        return;
    }

    // rebuild the boundary-cell mask when the mesh has changed;
    // amrStencil rebuilds its own caches internally
    if (boundaryCells_.size() != mesh_.nCells() || mesh_.changing())
    {
        boundaryCells_.setSize(mesh_.nCells());
        boundaryCells_ = false;
        markBoundaryCells(1);
    }

    interfaceCell_ = false;

    interfaceLabels_.clear();

    forAll(alpha1_, celli)
    {
        if (surfCellTol_ < alpha1_[celli] && alpha1_[celli] < 1-surfCellTol_)
        {
            interfaceCell_[celli] = true;
            interfaceLabels_.append(celli);
        }
    }
    interfaceNormal_.setSize(interfaceLabels_.size());
    centre_ = dimensionedVector("centre", dimLength, vector::zero);
    normal_ = dimensionedVector("normal", dimArea, vector::zero);

    gradSurf(alpha1_);

    forAll(interfaceLabels_, i)
    {
        const label celli = interfaceLabels_[i];
        if (mag(interfaceNormal_[i]) == 0)
        {
            continue;
        }

        sIterPLIC_.vofCutCell
        (
            celli,
            alpha1_[celli],
            isoFaceTol_,
            100,
            interfaceNormal_[i]
        );

        if (sIterPLIC_.cellStatus() == 0)
        {
            normal_[celli] = sIterPLIC_.surfaceArea();
            centre_[celli] = sIterPLIC_.surfaceCentre();
            if (mag(normal_[celli]) == 0)
            {
                normal_[celli] = vector::zero;
                centre_[celli] = vector::zero;
            }
        }
        else
        {
            normal_[celli] = vector::zero;
            centre_[celli] = vector::zero;
        }
    }
}


void Foam::reconstruction::mycCartesianAMR::mapAlphaField() const
{
    // without this line, we seem to get a race condition
    mesh_.C();

    cutCellPLIC cutCell(mesh_);

    forAll(normal_, celli)
    {
        if (mag(normal_[celli]) != 0)
        {
            vector n = normal_[celli]/mag(normal_[celli]);
            scalar cutValue = (centre_[celli] - mesh_.C()[celli]) & (n);
            cutCell.calcSubCell
            (
                celli,
                cutValue,
                n
            );
            alpha1_[celli] = cutCell.VolumeOfFluid();
        }
    }
    alpha1_.correctBoundaryConditions();
    alpha1_.oldTime() = alpha1_;
    alpha1_.oldTime().correctBoundaryConditions();
}


// ************************************************************************* //
