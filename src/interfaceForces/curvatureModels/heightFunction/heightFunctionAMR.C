/*---------------------------------------------------------------------------*\
            Copyright (c) 2021-2026, Asim Onder
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
    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "heightFunctionAMR.H"
#include "addToRunTimeSelectionTable.H"

#include "alphaContactAngleTwoPhaseFvPatchScalarField.H"
#include "mathematicalConstants.H"
#include "surfaceInterpolate.H"
#include "fvcDiv.H"
#include "fvcGrad.H"
#include "reconstructionSchemes.H"
#include "OFstream.H"
#include "IFstream.H"
#include "unitConversion.H"
#include "emptyPolyPatch.H"
#include "cyclicPolyPatch.H"
#include "symmetryPolyPatch.H"
#include "processorPolyPatch.H"
#include "wedgePolyPatch.H"
#include "DynamicList.H"
#include "reconstructedDistanceFunction.H"
#include "rdfCurvature.H"

namespace Foam
{
    defineTypeNameAndDebug(heightFunctionAMR, 0);
    addToRunTimeSelectionTable(curvatureModel, heightFunctionAMR, components);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::heightFunctionAMR::heightFunctionAMR
(
    const dictionary& dict,
    const volScalarField& alpha1,
    const volVectorField& U
)
:
    curvatureModel(typeName, dict, alpha1, U),
    mesh_(alpha1.mesh()),
    deltaN_("deltaN", 1e-8/pow(average(alpha1.mesh().V()), 1.0/3.0)),
    isEmpty_(false, false, false),
    emptyDirThickness_(0),
    nMax_(dict.subDict("heightFunctionParams").lookupOrDefault<label>("nMax", 3)),
    is2D_(false),
    interfaceTol_(1e-12),
    boundaryCells_(mesh_.nCells(), false),
    fallbackMethod_
    (
        dict.subDict("heightFunctionParams")
            .lookupOrDefault<word>("fallbackMethod", "rdf")
    )
{
    Info << "Curvature model: heightFunctionAMR." << endl;
    Info << "heightFunctionAMR: fallbackMethod = " << fallbackMethod_ << endl;

    const polyBoundaryMesh& bMesh = mesh_.boundaryMesh();
    forAll(bMesh, patchi)
    {
        const polyPatch& pp = bMesh[patchi];
        if (isA<emptyPolyPatch>(pp) || isA<wedgePolyPatch>(pp))
        {
            if (pp.size() > 0)
            {
                vector nf = pp.faceAreas()[0] / mag(pp.faceAreas()[0]);
                if (mag(nf.x()) > 0.5) isEmpty_.x() = true;
                else if (mag(nf.y()) > 0.5) isEmpty_.y() = true;
                else if (mag(nf.z()) > 0.5) isEmpty_.z() = true;
            }
        }
    }
    is2D_ = (isEmpty_.x() || isEmpty_.y() || isEmpty_.z());

    boundBox bb = mesh_.bounds();
    reduce(bb.min(), minOp<point>());
    reduce(bb.max(), maxOp<point>());
    if (isEmpty_.x()) emptyDirThickness_ = bb.max().x() - bb.min().x();
    else if (isEmpty_.y()) emptyDirThickness_ = bb.max().y() - bb.min().y();
    else if (isEmpty_.z()) emptyDirThickness_ = bb.max().z() - bb.min().z();

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
                vector nf = fAreas[0] / mag(fAreas[0]);
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

    scalar Lx = bb.max().x()-bb.min().x();
    scalar Ly = bb.max().y()-bb.min().y();
    scalar Lz = bb.max().z()-bb.min().z();

    Info << "heightFunctionAMR: isCyclic=("
         << isCyclic.x() << " " << isCyclic.y() << " " << isCyclic.z()
         << ") L=(" << Lx << " " << Ly << " " << Lz << ")" << endl;

    markBoundaryCells(nMax_);
    isCyclic_ = isCyclic;
    Lx_ = Lx; Ly_ = Ly; Lz_ = Lz;

    scalar dx = isEmpty_.x() ? emptyDirThickness_ : 0;
    scalar dy = isEmpty_.y() ? emptyDirThickness_ : 0;
    scalar dz = isEmpty_.z() ? emptyDirThickness_ : 0;
    amrStencilPtr_.reset
    (
        new amrStencil(mesh_, nMax_, nMax_, nMax_,
            dx, dy, dz, isEmpty_, isCyclic, Lx, Ly, Lz)
    );
}


void Foam::heightFunctionAMR::markBoundaryCells(label nLayers)
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


void Foam::heightFunctionAMR::correctContactAngle
(
    surfaceVectorField::Boundary& nHatb,
    const surfaceVectorField::Boundary& gradAlphaf
) const
{
    const fvMesh& mesh = alpha1_.mesh();
    const volScalarField::Boundary& abf = alpha1_.boundaryField();
    const fvBoundaryMesh& boundary = mesh.boundary();
    forAll(boundary, patchi)
    {
        if (isA<alphaContactAngleTwoPhaseFvPatchScalarField>(abf[patchi]))
        {
            alphaContactAngleTwoPhaseFvPatchScalarField& acap =
                const_cast<alphaContactAngleTwoPhaseFvPatchScalarField&>
                (refCast<const alphaContactAngleTwoPhaseFvPatchScalarField>(abf[patchi]));
            fvsPatchVectorField& nHatp = nHatb[patchi];
            const scalarField theta(degToRad()*acap.theta(U_.boundaryField()[patchi], nHatp));
            const vectorField nf(boundary[patchi].nf());
            const scalarField a12(nHatp & nf);
            const scalarField b1(cos(theta));
            scalarField b2(nHatp.size());
            forAll(b2, facei)
                b2[facei] = cos(acos(a12[facei]) - theta[facei]);
            const scalarField det(1.0 - a12*a12);
            scalarField a((b1-a12*b2)/det);
            scalarField b((b2-a12*b1)/det);
            nHatp = a*nf + b*nHatp;
            nHatp /= (mag(nHatp) + deltaN_.value());
            acap.gradient() = (nf & nHatp)*mag(gradAlphaf[patchi]);
            acap.evaluate();
        }
    }
}


void Foam::heightFunctionAMR::calculateK()
{
    const fvMesh& mesh = alpha1_.mesh();

    Map<scalar> alphaIJK;

    forAll(K_, celli)
        K_[celli] = 0;

    reconstructionSchemes& surf =
        mesh.lookupObjectRef<reconstructionSchemes>("reconstructionScheme");

    const boolList& interfaceCells = surf.interfaceCell();
    const volVectorField& faceNormal = surf.normal();

    // --- T1: getZoneField (comm + hash build) ---
    scalar t1a = mesh.time().elapsedCpuTime();
    amrStencilPtr_->getZoneField(interfaceCells, alphaIJK, alpha1_);
    scalar t1b = mesh.time().elapsedCpuTime();

    const Map<scalar>& aHash = amrStencilPtr_->alphaHash();
    
    

    // --- Self-test: verify hash correctness at fine-level seed positions ---
    {
        scalar testFineV = GREAT;
        forAll(interfaceCells, ci)
        {
            if (interfaceCells[ci])
                testFineV = min(testFineV, mesh_.V()[ci]);
        }
        reduce(testFineV, minOp<scalar>());

        label nSeedMiss = 0;
        label nSeedWrong = 0;
        label nSeedOK = 0;
        forAll(interfaceCells, celli)
        {
            if (!interfaceCells[celli]) continue;
            if (mesh_.V()[celli] > 1.5 * testFineV) continue;  // skip coarse
            point p = amrStencilPtr_->cyclicWrap(mesh.C()[celli]);
            const label key = amrStencilPtr_->posKey(p);
            const auto iter = aHash.cfind(key);
            if (!iter.good())
                nSeedMiss++;
            else if (mag(iter.val() - alpha1_[celli]) > 1e-10)
                nSeedWrong++;
            else
                nSeedOK++;
        }
        reduce(nSeedMiss, sumOp<label>());
        reduce(nSeedWrong, sumOp<label>());
        reduce(nSeedOK, sumOp<label>());
        Info<< "heightFunctionAMR: hashTest OK=" << nSeedOK
            << " miss=" << nSeedMiss
            << " wrong=" << nSeedWrong << endl;
    }

    label nIncomplete = 0;
    label nBoundaryFallback = 0;
    label nCoarseLevel = 0;
    label nHFCells = 0;
    DynamicList<label> fallbackCells;

    // Compute finest interface cell volume for seed-level check
    scalar fineV = GREAT;
    forAll(interfaceCells, ci)
    {
        if (interfaceCells[ci])
            fineV = min(fineV, mesh_.V()[ci]);
    }
    reduce(fineV, minOp<scalar>());
    const scalar seedVolTol = 1.5 * fineV;

    // --- T2: direction-aware HF via direct hash lookup ---
    scalar t2a = mesh.time().elapsedCpuTime();

    forAll(interfaceCells, celli)
    {
        vector n = faceNormal[celli];
        if (!interfaceCells[celli] || mag(n) == 0) continue;

        nHFCells++;

        label dir = -1;
        if (mag(n.x()) > 0.0 && !isEmpty_.x()) dir = 0;
        if (mag(n.y()) > mag(n.x()) && !isEmpty_.y()) dir = 1;
        if (mag(n.z()) > mag(n.y()) && mag(n.z()) > mag(n.x())
            && !isEmpty_.z()) dir = 2;
        if (dir == -1) continue;

        if (celli < boundaryCells_.size() && boundaryCells_[celli])
        {
            fallbackCells.append(celli);
            nBoundaryFallback++;
            continue;
        }

        // Skip coarse-level interface cells (hash is at finest level)
        if (mesh_.V()[celli] > seedVolTol)
        {
            fallbackCells.append(celli);
            nCoarseLevel++;
            continue;
        }

        const scalar V = mesh_.V()[celli];
        scalar dl;
        if (is2D_)
            dl = Foam::sqrt(V / emptyDirThickness_);
        else
            dl = Foam::cbrt(V);

        const point seedCC = mesh.C()[celli];
        bool incomplete = false;

        if (is2D_)
        {
            List<scalar> H(3, 0.0);

            label sweepAxis = dir;
            label crossAxis = -1;

            if (isEmpty_.x())
            {
                if (dir == 1) { sweepAxis = 1; crossAxis = 2; }
                else          { sweepAxis = 2; crossAxis = 1; }
            }
            else if (isEmpty_.y())
            {
                if (dir == 0) { sweepAxis = 0; crossAxis = 2; }
                else          { sweepAxis = 2; crossAxis = 0; }
            }
            else
            {
                if (dir == 0) { sweepAxis = 0; crossAxis = 1; }
                else          { sweepAxis = 1; crossAxis = 0; }
            }

            for (int c = -1; c < 2 && !incomplete; c++)
            for (int s = -nMax_; s <= nMax_ && !incomplete; s++)
            {
                point eCC = seedCC;
                eCC[sweepAxis] += s * dl;
                eCC[crossAxis] += c * dl;
                eCC = amrStencilPtr_->cyclicWrap(eCC);

                const label key = amrStencilPtr_->posKey(eCC);
                const auto iter = aHash.cfind(key);

                if (iter.good())
                    H[c+1] += iter.val();
                else
                    incomplete = true;
            }

            if (incomplete)
            {
                nIncomplete++;
                fallbackCells.append(celli);
                continue;
            }

            H = -H * dl;

            scalar Ht = (H[2]-H[0]) / 2.0 / dl;
            scalar Htt = (H[2]-2.0*H[1]+H[0]) / dl / dl;
            K_[celli] = Htt / Foam::pow(1.0+Ht*Ht, 1.5);
        }
        else
        {
            List<scalar> H(9, 0.0);

            label a1, a2c;
            if (dir == 0)      { a1 = 1; a2c = 2; }
            else if (dir == 1) { a1 = 0; a2c = 2; }
            else               { a1 = 0; a2c = 1; }

            for (int c1 = -1; c1 < 2 && !incomplete; c1++)
            for (int c2 = -1; c2 < 2 && !incomplete; c2++)
            for (int s = -nMax_; s <= nMax_ && !incomplete; s++)
            {
                point eCC = seedCC;
                eCC[dir] += s * dl;
                eCC[a1]  += c1 * dl;
                eCC[a2c] += c2 * dl;
                eCC = amrStencilPtr_->cyclicWrap(eCC);

                const label key = amrStencilPtr_->posKey(eCC);
                const auto iter = aHash.cfind(key);

                if (iter.good())
                    H[c1+1 + 3*(c2+1)] += iter.val();
                else
                    incomplete = true;
            }

            if (incomplete)
            {
                nIncomplete++;
                fallbackCells.append(celli);
                continue;
            }

            H = -H * dl;

            double H2[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
            for (int i = -1; i < 2; i++)
                for (int j = -1; j < 2; j++)
                    H2[i+1][j+1] = H[i+1 + 3*(j+1)];

            scalar Ht = (H2[2][1]-H2[0][1]) / 2.0 / dl;
            scalar Hb = (H2[1][2]-H2[1][0]) / 2.0 / dl;
            scalar Htt = (H2[2][1]-2.0*H2[1][1]+H2[0][1]) / dl / dl;
            scalar Hbb = (H2[1][2]-2.0*H2[1][1]+H2[1][0]) / dl / dl;
            scalar Htb = (H2[2][2]-H2[2][0]-H2[0][2]+H2[0][0]) / 4.0 / dl / dl;
            K_[celli] =
                (Htt+Hbb+Htt*Hb*Hb+Hbb*Ht*Ht-2.0*Htb*Ht*Hb)
                / Foam::pow(1.0+Ht*Ht+Hb*Hb, 1.5);
        }
    }

    scalar t2b = mesh.time().elapsedCpuTime();

    // --- T3: fallback ---
    scalar t3a = mesh.time().elapsedCpuTime();

    label nFallbackGlobal = fallbackCells.size();
    reduce(nFallbackGlobal, sumOp<label>());

    if (nFallbackGlobal > 0)
    {
        if (fallbackMethod_ == "rdf")
        {
            tmp<volScalarField> tKRDF =
                rdfCurvature::computeKRDF(mesh, alpha1_);
            const volScalarField& KRDF = tKRDF();
            forAll(fallbackCells, i)
                K_[fallbackCells[i]] = KRDF[fallbackCells[i]];
        }
        else if (fallbackMethod_ == "divGradAlpha")
        {
            volVectorField gradAlpha(fvc::grad(alpha1_, "nHat"));
            gradAlpha /= (mag(gradAlpha)
                + dimensionedScalar("small", dimless, 1e-14));
            gradAlpha.correctBoundaryConditions();
            const volScalarField KDivGrad(-tr(fvc::grad(gradAlpha, "nHat")));

            forAll(fallbackCells, i)
                K_[fallbackCells[i]] = KDivGrad[fallbackCells[i]];
        }
        else if (fallbackMethod_ == "interpolate")
        {
            // Average K from face-neighboring cells that have valid K
            const labelUList& own = mesh.faceOwner();
            const labelUList& nei = mesh.faceNeighbour();

            // Mark fallback cells
            boolList isFallback(mesh.nCells(), false);
            forAll(fallbackCells, i)
                isFallback[fallbackCells[i]] = true;

            forAll(fallbackCells, i)
            {
                const label celli = fallbackCells[i];
                const cell& c = mesh.cells()[celli];
                scalar sumK = 0;
                label nK = 0;

                forAll(c, fi)
                {
                    const label faceI = c[fi];
                    if (faceI >= mesh.nInternalFaces()) continue;

                    const label nbr = (own[faceI] == celli)
                        ? nei[faceI] : own[faceI];

                    if (!isFallback[nbr] && K_[nbr] != 0)
                    {
                        sumK += K_[nbr];
                        nK++;
                    }
                }

                if (nK > 0)
                    K_[celli] = sumK / nK;
                // else: stays 0
            }
        }
        // else: fallbackMethod_ == "none" -> K stays 0
    }

    scalar t3b = mesh.time().elapsedCpuTime();

    // --- T4: diagnostics + targeted face interpolation ---
    scalar t4a = mesh.time().elapsedCpuTime();

    reduce(nHFCells, sumOp<label>());
    reduce(nIncomplete, sumOp<label>());
    reduce(nBoundaryFallback, sumOp<label>());
    reduce(nCoarseLevel, sumOp<label>());

    const label nFallback = nIncomplete + nBoundaryFallback + nCoarseLevel;
    if (nFallback > 0)
    {
        Info<< "heightFunctionAMR: " << nHFCells << " interface cells, "
            << nFallback << " fell back (" << fallbackMethod_ << ")"
            << " (boundary=" << nBoundaryFallback
            << " incomplete=" << nIncomplete
            << " coarseLevel=" << nCoarseLevel << ")" << endl;
    }

    // Zero K at near-pure interface cells
    forAll(interfaceCells, iCell)
    {
        if (interfaceCells[iCell])
        {
            if (alpha1_[iCell] < interfaceTol_
                || alpha1_[iCell] > 1 - interfaceTol_)
                K_[iCell] = 0.0;
        }
    }

    // --- Targeted Kf with processor boundary support ---
    K_.correctBoundaryConditions();
    Kf_ = dimensionedScalar("zero", Kf_.dimensions(), 0);

    {
        const labelUList& own = mesh.faceOwner();
        const labelUList& nei = mesh.faceNeighbour();
        const scalarField& weights = mesh.weights();

        // Internal faces adjacent to interface cells
        forAll(interfaceCells, celli)
        {
            if (!interfaceCells[celli]) continue;
            const cell& c = mesh.cells()[celli];
            forAll(c, fi)
            {
                const label faceI = c[fi];
                if (faceI >= mesh.nInternalFaces()) continue;
                const scalar K1 = K_[own[faceI]];
                const scalar K2 = K_[nei[faceI]];
                if (K1 != 0 && K2 != 0)
                    Kf_[faceI] = weights[faceI]*K1
                        + (1.0 - weights[faceI])*K2;
                else if (K1 != 0)
                    Kf_[faceI] = K1;
                else if (K2 != 0)
                    Kf_[faceI] = K2;
            }
        }

        // Boundary faces: use K_.boundaryField() for neighbor values
        const scalar dl_ = is2D_
            ? Foam::sqrt(average(mesh.V()).value() / emptyDirThickness_)
            : Foam::cbrt(average(mesh.V()).value());
        const scalar thr = 1e-14 / dl_;

        forAll(mesh.boundaryMesh(), iPatch)
        {
            const word pType = Kf_.boundaryField()[iPatch].type();
            if (pType == "empty") continue;

            const polyPatch& cPatch = mesh.boundaryMesh()[iPatch];
            const labelUList& faceCells = cPatch.faceCells();
            const scalarField& Kb = K_.boundaryField()[iPatch];
            const scalarField& wt = mesh.weights().boundaryField()[iPatch];

            forAll(faceCells, iFace)
            {
                const label iCell = faceCells[iFace];
                const scalar Klocal = K_[iCell];
                const scalar Kremote = Kb[iFace];

                if (Klocal != 0 && Kremote != 0)
                {
                    Kf_.boundaryFieldRef()[iPatch][iFace] =
                        wt[iFace]*Klocal + (1.0 - wt[iFace])*Kremote;
                }
                else if (Klocal != 0)
                {
                    Kf_.boundaryFieldRef()[iPatch][iFace] = Klocal;
                }
                else if (Kremote != 0)
                {
                    Kf_.boundaryFieldRef()[iPatch][iFace] = Kremote;
                }
            }
        }
    }

    scalar t4b = mesh.time().elapsedCpuTime();

    scalar dt1=t1b-t1a, dt2=t2b-t2a, dt3=t3b-t3a, dt4=t4b-t4a;
    reduce(dt1, maxOp<scalar>());
    reduce(dt2, maxOp<scalar>());
    reduce(dt3, maxOp<scalar>());
    reduce(dt4, maxOp<scalar>());

    Info << "heightFunctionAMR: T1(getZoneField)=" << dt1
         << "s T2(HFloop)=" << dt2
         << "s T3(fallback)=" << dt3
         << "s T4(faceInterp)=" << dt4 << "s" << endl;
}

// ************************************************************************* //
