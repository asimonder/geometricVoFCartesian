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
    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "rdfCurvature.H"
#include "addToRunTimeSelectionTable.H"

#include "alphaContactAngleTwoPhaseFvPatchScalarField.H"
#include "mathematicalConstants.H"
#include "surfaceInterpolate.H"
#include "fvcDiv.H"
#include "fvcGrad.H"
#include "unitConversion.H"
#include "reconstructionSchemes.H"
#include "zoneDistribute.H"
#include "processorCyclicPolyPatch.H"
#include "indexedOctree.H"
#include "treeDataPoint.H"
#include "Random.H"
#include "reconstructedDistanceFunction.H"

namespace Foam
{
    defineTypeNameAndDebug(rdfCurvature, 0);
    addToRunTimeSelectionTable(curvatureModel, rdfCurvature, components);
}


Foam::rdfCurvature::rdfCurvature
(
    const dictionary& dict,
    const volScalarField& alpha1,
    const volVectorField& U
)
:
    curvatureModel(typeName, dict, alpha1, U),
    deltaN_
    (
        "deltaN",
        1e-8/pow(average(alpha1.mesh().V()), 1.0/3.0)
    ),
    nBand_
    (
        dict.subOrEmptyDict("rdfCurvatureParams")
            .lookupOrDefault<label>("nBand", 4)
     ),
    interfaceTol_(1e-12)
{
    Info << "Curvature model: rdfCurvature." << endl;
}


void Foam::rdfCurvature::correctContactAngle
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
                (
                    refCast<const alphaContactAngleTwoPhaseFvPatchScalarField>
                    (
                        abf[patchi]
                    )
                );

            fvsPatchVectorField& nHatp = nHatb[patchi];
            const scalarField theta
            (
                degToRad() * acap.theta(U_.boundaryField()[patchi], nHatp)
            );

            const vectorField nf(boundary[patchi].nf());

            const scalarField a12(nHatp & nf);
            const scalarField b1(cos(theta));

            scalarField b2(nHatp.size());
            forAll(b2, facei)
            {
                b2[facei] = cos(acos(a12[facei]) - theta[facei]);
            }

            const scalarField det(1.0 - a12*a12);

            scalarField a((b1 - a12*b2)/det);
            scalarField b((b2 - a12*b1)/det);

            nHatp = a*nf + b*nHatp;
            nHatp /= (mag(nHatp) + deltaN_.value());

            acap.gradient() = (nf & nHatp)*mag(gradAlphaf[patchi]);
            acap.evaluate();
        }
    }
}


Foam::tmp<Foam::volScalarField> Foam::rdfCurvature::computeKRDF
(
    const fvMesh& mesh,
    const volScalarField& alpha1
)
{
    reconstructionSchemes& surf =
        mesh.lookupObjectRef<reconstructionSchemes>("reconstructionScheme");

    surf.reconstruct(false);

    const volVectorField& faceCentre = surf.centre();
    const volVectorField& faceNormal = surf.normal();
    const boolList& interfaceCellsRef = surf.interfaceCell();

    boolList interfaceCells(interfaceCellsRef.size(), false);
    forAll(interfaceCellsRef, ci)
    {
        if (interfaceCellsRef[ci] && mag(faceNormal[ci]) != 0)
        {
            interfaceCells[ci] = true;
        }
    }

    // Band
    reconstructedDistanceFunction RDF(mesh);
    RDF.markCellsNearSurf(interfaceCells, 4);
    const boolList& band = RDF.nextToInterface();

    // Collect local interface centres/normals
    DynamicField<vector> cLocal(mesh.nCells());
    DynamicField<vector> nLocal(mesh.nCells());
    forAll(interfaceCells, ci)
    {
        if (interfaceCells[ci])
        {
            cLocal.append(faceCentre[ci]);
            nLocal.append(faceNormal[ci]);
        }
    }

    // Neighbor-only exchange
    DynamicField<vector> cAll(cLocal);
    DynamicField<vector> nAll(nLocal);

    if (Pstream::parRun())
    {
        HashSet<label> nbrSet;
        {
            const polyBoundaryMesh& patches = mesh.boundaryMesh();
            forAll(patches, patchi)
            {
                const polyPatch& pp = patches[patchi];
                if (isA<processorPolyPatch>(pp))
                {
                    const processorPolyPatch& procPatch =
                        refCast<const processorPolyPatch>(pp);
                    nbrSet.insert(procPatch.neighbProcNo());
                }
            }
        }
        const labelList partners(nbrSet.toc());

        // Round 1
        vectorField r1C, r1N;
        labelList r1Src;
        {
            PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);
            for (const label p : partners)
            {
                UOPstream toNbr(p, pBufs);
                toNbr << vectorField(cLocal) << vectorField(nLocal);
            }
            pBufs.finishedSends();

            DynamicField<vector> tmpC, tmpN;
            DynamicList<label> tmpS;
            for (const label p : partners)
            {
                vectorField recvC, recvN;
                UIPstream fromNbr(p, pBufs);
                fromNbr >> recvC >> recvN;
                forAll(recvC, i)
                {
                    tmpC.append(recvC[i]);
                    tmpN.append(recvN[i]);
                    tmpS.append(p);
                }
            }
            r1C = tmpC;
            r1N = tmpN;
            r1Src = tmpS;
        }
        cAll.append(r1C);
        nAll.append(r1N);

        // Round 2: relay for 2-hop
        {
            PstreamBuffers pBufs2(Pstream::commsTypes::nonBlocking);
            for (const label p : partners)
            {
                DynamicField<vector> sendC, sendN;
                forAll(r1C, ri)
                {
                    if (r1Src[ri] != p)
                    {
                        sendC.append(r1C[ri]);
                        sendN.append(r1N[ri]);
                    }
                }
                UOPstream toNbr(p, pBufs2);
                toNbr << vectorField(sendC) << vectorField(sendN);
            }
            pBufs2.finishedSends();
            for (const label p : partners)
            {
                vectorField recvC, recvN;
                UIPstream fromNbr(p, pBufs2);
                fromNbr >> recvC >> recvN;
                cAll.append(recvC);
                nAll.append(recvN);
            }
        }
    }

    // Cyclic images: GLOBAL domain bounds
    boundBox globalBB = mesh.bounds();
    reduce(globalBB.min(), minOp<point>());
    reduce(globalBB.max(), maxOp<point>());
    const vector L = globalBB.max() - globalBB.min();

    Vector<bool> isCyclic(false, false, false);
    {
        const polyBoundaryMesh& patches = mesh.boundaryMesh();
        forAll(patches, patchi)
        {
            const polyPatch& pp = patches[patchi];
            if (isA<cyclicPolyPatch>(pp) || isA<processorCyclicPolyPatch>(pp))
            {
                if (pp.size() == 0) continue;
                const vector nrm = pp.faceNormals()[0];
                if (mag(nrm.x()) > 0.9) isCyclic.x() = true;
                if (mag(nrm.y()) > 0.9) isCyclic.y() = true;
                if (mag(nrm.z()) > 0.9) isCyclic.z() = true;
            }
        }
        reduce(isCyclic.x(), orOp<bool>());
        reduce(isCyclic.y(), orOp<bool>());
        reduce(isCyclic.z(), orOp<bool>());
    }

    DynamicField<vector> cPool(cAll);
    DynamicField<vector> nPool(nAll);

    for (label sx = -1; sx <= 1; sx++)
    {
        if (sx != 0 && !isCyclic.x()) continue;
        for (label sy = -1; sy <= 1; sy++)
        {
            if (sy != 0 && !isCyclic.y()) continue;
            for (label sz = -1; sz <= 1; sz++)
            {
                if (sz != 0 && !isCyclic.z()) continue;
                if (sx == 0 && sy == 0 && sz == 0) continue;
                const vector shift(sx*L.x(), sy*L.y(), sz*L.z());
                forAll(cAll, k)
                {
                    cPool.append(cAll[k] + shift);
                    nPool.append(nAll[k]);
                }
            }
        }
    }

    // Octree over local + neighbor data (NOT global)
    Random rndGen(1234567);
    treeBoundBox bb(treeBoundBox(cPool).extend(rndGen, 1e-4));
    bb.min() -= point(1e-8, 1e-8, 1e-8);
    bb.max() += point(1e-8, 1e-8, 1e-8);

    indexedOctree<treeDataPoint> surfaceTree
    (
        treeDataPoint(cPool),
        bb, 8, 10, 3.0
    );

    // RDF via nearest-plane projection (single code path)
    volScalarField& RDFvol = RDF;
    forAll(RDFvol, ci)
    {
        if (!band[ci])
        {
            RDFvol[ci] = 0;
            continue;
        }

        const point p = mesh.C()[ci];
        pointIndexHit pHit = surfaceTree.findNearest(p, GREAT);
        const label iMin = pHit.index();

        if (iMin >= 0)
        {
            vector nHat = -nPool[iMin] / mag(nPool[iMin]);
            RDFvol[ci] = (cPool[iMin] - p) & nHat;
        }
        else
        {
            RDFvol[ci] = 0;
        }
    }
    RDFvol.correctBoundaryConditions();

    // Curvature
    volVectorField gradRDF(fvc::grad(RDFvol, "nHat"));
    dimensionedScalar deltaN("deltaN", dimless, 1e-14);
    gradRDF /= (mag(gradRDF) + deltaN);
    gradRDF.correctBoundaryConditions();

    tmp<volScalarField> tKRDF
    (
        new volScalarField
        (
            IOobject
            (
                "KRDF",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            -tr(fvc::grad(gradRDF, "nHat"))
        )
    );

    return tKRDF;
}

void Foam::rdfCurvature::calculateK()
{
    const fvMesh& mesh = alpha1_.mesh();

    reconstructedDistanceFunction RDF(mesh);
    
    // --- 1. Reconstruction data ---
    reconstructionSchemes& surf =
        mesh.lookupObjectRef<reconstructionSchemes>("reconstructionScheme");

    surf.reconstruct(false);

    const volVectorField& faceCentre = surf.centre();
    const volVectorField& faceNormal = surf.normal();
    const boolList& interfaceCellsRef = surf.interfaceCell();

    boolList interfaceCells(interfaceCellsRef.size(), false);
    forAll(interfaceCellsRef, ci)
    {
        if (interfaceCellsRef[ci] && mag(faceNormal[ci]) != 0)
        {
            interfaceCells[ci] = true;
        }
    }

    // --- 2. Band ---
    RDF.markCellsNearSurf(interfaceCells, nBand_);
    const boolList& band = RDF.nextToInterface();

    // --- 3. Collect local interface centres and normals ---
    DynamicField<vector> cLocal(mesh.nCells());
    DynamicField<vector> nLocal(mesh.nCells());
    forAll(interfaceCells, ci)
    {
        if (interfaceCells[ci])
        {
            cLocal.append(faceCentre[ci]);
            nLocal.append(faceNormal[ci]);
        }
    }

    // Face-adjacent partners (reused for K exchange)
    labelList partners;
    if (Pstream::parRun())
    {
        HashSet<label> nbrSet;
        const polyBoundaryMesh& bMesh = mesh.boundaryMesh();
        forAll(bMesh, patchi)
        {
            const polyPatch& pp = bMesh[patchi];
            if (isA<processorPolyPatch>(pp))
            {
                const processorPolyPatch& procPatch =
                    refCast<const processorPolyPatch>(pp);
                nbrSet.insert(procPatch.neighbProcNo());
            }
        }
        partners = nbrSet.toc();
    }

    // --- 3b. Neighbor-only exchange of centres/normals ---
    DynamicField<vector> cAll(cLocal);
    DynamicField<vector> nAll(nLocal);

    if (Pstream::parRun())
    {
        // Round 1
        vectorField r1C, r1N;
        labelList r1Src;
        {
            PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);
            for (const label p : partners)
            {
                UOPstream toNbr(p, pBufs);
                toNbr << vectorField(cLocal) << vectorField(nLocal);
            }
            pBufs.finishedSends();

            DynamicField<vector> tmpC, tmpN;
            DynamicList<label> tmpS;
            for (const label p : partners)
            {
                vectorField recvC, recvN;
                UIPstream fromNbr(p, pBufs);
                fromNbr >> recvC >> recvN;
                forAll(recvC, i)
                {
                    tmpC.append(recvC[i]);
                    tmpN.append(recvN[i]);
                    tmpS.append(p);
                }
            }
            r1C = tmpC;
            r1N = tmpN;
            r1Src = tmpS;
        }
        cAll.append(r1C);
        nAll.append(r1N);

        // Round 2: relay for 2-hop
        {
            PstreamBuffers pBufs2(Pstream::commsTypes::nonBlocking);
            for (const label p : partners)
            {
                DynamicField<vector> sendC, sendN;
                forAll(r1C, ri)
                {
                    if (r1Src[ri] != p)
                    {
                        sendC.append(r1C[ri]);
                        sendN.append(r1N[ri]);
                    }
                }
                UOPstream toNbr(p, pBufs2);
                toNbr << vectorField(sendC) << vectorField(sendN);
            }
            pBufs2.finishedSends();
            for (const label p : partners)
            {
                vectorField recvC, recvN;
                UIPstream fromNbr(p, pBufs2);
                fromNbr >> recvC >> recvN;
                cAll.append(recvC);
                nAll.append(recvN);
            }
        }
    }

    // --- 4. Cyclic image copies (GLOBAL domain bounds) ---
    boundBox globalBB = mesh.bounds();
    reduce(globalBB.min(), minOp<point>());
    reduce(globalBB.max(), maxOp<point>());
    const vector L = globalBB.max() - globalBB.min();

    Vector<bool> isCyclic(false, false, false);
    {
        const polyBoundaryMesh& patches = mesh.boundaryMesh();
        forAll(patches, patchi)
        {
            const polyPatch& pp = patches[patchi];
            if (isA<cyclicPolyPatch>(pp) || isA<processorCyclicPolyPatch>(pp))
            {
                if (pp.size() == 0) continue;
                const vector nrm = pp.faceNormals()[0];
                if (mag(nrm.x()) > 0.9) isCyclic.x() = true;
                if (mag(nrm.y()) > 0.9) isCyclic.y() = true;
                if (mag(nrm.z()) > 0.9) isCyclic.z() = true;
            }
        }
        reduce(isCyclic.x(), orOp<bool>());
        reduce(isCyclic.y(), orOp<bool>());
        reduce(isCyclic.z(), orOp<bool>());
    }

    DynamicField<vector> cPool(cAll);
    DynamicField<vector> nPool(nAll);

    for (label sx = -1; sx <= 1; sx++)
    {
        if (sx != 0 && !isCyclic.x()) continue;
        for (label sy = -1; sy <= 1; sy++)
        {
            if (sy != 0 && !isCyclic.y()) continue;
            for (label sz = -1; sz <= 1; sz++)
            {
                if (sz != 0 && !isCyclic.z()) continue;
                if (sx == 0 && sy == 0 && sz == 0) continue;

                const vector shift(sx*L.x(), sy*L.y(), sz*L.z());
                forAll(cAll, k)
                {
                    cPool.append(cAll[k] + shift);
                    nPool.append(nAll[k]);
                }
            }
        }
    }

    // --- 5. Octree over local + neighbor data ---
    Random rndGen(1234567);
    treeBoundBox bb(treeBoundBox(cPool).extend(rndGen, 1e-4));
    bb.min() -= point(1e-8, 1e-8, 1e-8);
    bb.max() += point(1e-8, 1e-8, 1e-8);

    indexedOctree<treeDataPoint> surfaceTree
    (
        treeDataPoint(cPool),
        bb,
        8,
        10,
        3.0
    );

    // --- 6. RDF via nearest-plane projection ---
    volScalarField& RDFvol = RDF;
    forAll(RDFvol, ci)
    {
        if (!band[ci])
        {
            RDFvol[ci] = 0;
            continue;
        }

        const point p = mesh.C()[ci];
        pointIndexHit pHit = surfaceTree.findNearest(p, GREAT);
        const label iMin = pHit.index();

        if (iMin >= 0)
        {
            vector nHat = -nPool[iMin] / mag(nPool[iMin]);
            RDFvol[ci] = (cPool[iMin] - p) & nHat;
        }
        else
        {
            RDFvol[ci] = 0;
        }
    }
    RDFvol.correctBoundaryConditions();

    // --- 7. Curvature from trace of normalized gradRDF ---
    volVectorField gradRDF(fvc::grad(RDFvol, "nHat"));
    dimensionedScalar deltaN("deltaN", dimless, 1e-14);
    gradRDF /= (mag(gradRDF) + deltaN);
    gradRDF.correctBoundaryConditions();

    const volScalarField KRDF(-tr(fvc::grad(gradRDF, "nHat")));

    // --- 8. K at interface cells only ---
    forAll(K_, ci) K_[ci] = 0;
    forAll(interfaceCells, ci)
    {
        if (interfaceCells[ci])
        {
            K_[ci] = KRDF[ci];
        }
    }

    // --- 9. Extend K to band cells from nearest interface cell ---
    DynamicField<scalar> kLocal;
    forAll(interfaceCells, ci)
    {
        if (interfaceCells[ci])
        {
            kLocal.append(K_[ci]);
        }
    }

    // Neighbor-only exchange of K values (same partners, same order)
    DynamicField<scalar> kAll(kLocal);

    if (Pstream::parRun())
    {
        // Round 1
        scalarField r1K;
        labelList r1KSrc;
        {
            PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);
            for (const label p : partners)
            {
                UOPstream toNbr(p, pBufs);
                toNbr << scalarField(kLocal);
            }
            pBufs.finishedSends();

            DynamicField<scalar> tmpK;
            DynamicList<label> tmpS;
            for (const label p : partners)
            {
                scalarField recvK;
                UIPstream fromNbr(p, pBufs);
                fromNbr >> recvK;
                forAll(recvK, i)
                {
                    tmpK.append(recvK[i]);
                    tmpS.append(p);
                }
            }
            r1K = tmpK;
            r1KSrc = tmpS;
        }
        kAll.append(r1K);

        // Round 2
        {
            PstreamBuffers pBufs2(Pstream::commsTypes::nonBlocking);
            for (const label p : partners)
            {
                DynamicField<scalar> sendK;
                forAll(r1K, ri)
                {
                    if (r1KSrc[ri] != p)
                    {
                        sendK.append(r1K[ri]);
                    }
                }
                UOPstream toNbr(p, pBufs2);
                toNbr << scalarField(sendK);
            }
            pBufs2.finishedSends();
            for (const label p : partners)
            {
                scalarField recvK;
                UIPstream fromNbr(p, pBufs2);
                fromNbr >> recvK;
                kAll.append(recvK);
            }
        }
    }

    DynamicField<scalar> kPool(kAll);
    for (label sx = -1; sx <= 1; sx++)
    {
        if (sx != 0 && !isCyclic.x()) continue;
        for (label sy = -1; sy <= 1; sy++)
        {
            if (sy != 0 && !isCyclic.y()) continue;
            for (label sz = -1; sz <= 1; sz++)
            {
                if (sz != 0 && !isCyclic.z()) continue;
                if (sx == 0 && sy == 0 && sz == 0) continue;
                forAll(kAll, k) kPool.append(kAll[k]);
            }
        }
    }

    // Reuse octree for K extension (same cPool geometry)
    forAll(K_, ci)
    {
        if (!band[ci] || interfaceCells[ci]) continue;

        const point p = mesh.C()[ci];
        pointIndexHit pHit = surfaceTree.findNearest(p, GREAT);
        const label iMin = pHit.index();

        if (iMin >= 0) K_[ci] = kPool[iMin];
    }

    K_.correctBoundaryConditions();
    Kf_ = fvc::interpolate(K_);
}
