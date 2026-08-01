/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2025 Asim Onder
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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

#include "amrStencil.H"
#include "PstreamReduceOps.H"
#include <queue>

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::amrStencil::amrStencil
(
    const fvMesh& mesh,
    int iMax,
    int jMax,
    int kMax,
    scalar dx,
    scalar dy,
    scalar dz,
    const Vector<bool>& isEmpty,
    const Vector<bool>& isCyclic,
    scalar Lx,
    scalar Ly,
    scalar Lz
)
:
    uniformStencil(SubclassTag{}, mesh, iMax, jMax, kMax, isEmpty, dx),
    mesh_(mesh),
    C_(mesh.C()),
    globalNumbering_(mesh.nCells() + mesh.nBoundaryFaces()),
    topoVersion_(-1),
    nCellsAtUpdate_(0),
    isCyclic_(isCyclic),
    Lx_(Lx),
    Ly_(Ly),
    Lz_(Lz),
    procMapsBuilt_(false)
{
    // Build neighbour proc list from processor patches
    if (Pstream::parRun())
    {
        const polyBoundaryMesh& bMesh = mesh_.boundaryMesh();
        HashSet<label> procNeighbours;

        forAll(bMesh, patchi)
        {
            const polyPatch& pp = bMesh[patchi];
            if (isA<processorPolyPatch>(pp))
            {
                const processorPolyPatch& procPatch =
                    refCast<const processorPolyPatch>(pp);
                procNeighbours.insert(procPatch.neighbProcNo());
            }
        }
        neighbourProcs_ = procNeighbours.toc();
    }

    // Build all caches and set topology version
    update();

    Info << "amrStencil: iMax=" << iMax_
         << " jMax=" << jMax_
         << " kMax=" << kMax_
         << " N=" << N_
         << " nNeighbourProcs=" << neighbourProcs_.size() << endl;
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::amrStencil::buildProcMaps()
{
    const polyBoundaryMesh& bMesh = mesh_.boundaryMesh();

    procPatchRemoteGId_.setSize(bMesh.size());
    procPatchRemoteCC_.setSize(bMesh.size());

    if (!Pstream::parRun())
    {
        procMapsBuilt_ = true;
        return;
    }

    PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

    forAll(bMesh, patchi)
    {
        const polyPatch& pp = bMesh[patchi];

        if (isA<processorPolyPatch>(pp))
        {
            const processorPolyPatch& procPatch =
                refCast<const processorPolyPatch>(pp);
            const labelUList& faceCells = procPatch.faceCells();

            labelList myGIds(faceCells.size());
            pointField myCCs(faceCells.size());

            forAll(faceCells, fi)
            {
                myGIds[fi] = globalNumbering_.toGlobal(faceCells[fi]);
                myCCs[fi] = C_[faceCells[fi]];
            }

            UOPstream toNbr(procPatch.neighbProcNo(), pBufs);
            toNbr << myGIds << myCCs;
        }
    }

    pBufs.finishedSends();

    forAll(bMesh, patchi)
    {
        const polyPatch& pp = bMesh[patchi];

        if (isA<processorPolyPatch>(pp))
        {
            const processorPolyPatch& procPatch =
                refCast<const processorPolyPatch>(pp);

            UIPstream fromNbr(procPatch.neighbProcNo(), pBufs);
            fromNbr >> procPatchRemoteGId_[patchi]
                    >> procPatchRemoteCC_[patchi];
        }
    }

    procMapsBuilt_ = true;
}


Foam::labelList Foam::amrStencil::buildStencilTopology
(
    const label seedCelli,
    List<labelHashSet>& commNeeded
) const
{
    // Result: globalId at each stencil position, -1 if not found
    labelList topo(N_, -1);

    const point seedCC = C_[seedCelli];
    const labelUList& own = mesh_.faceOwner();
    const labelUList& nei = mesh_.faceNeighbour();
    const label nInternalFaces = mesh_.nInternalFaces();
    const polyBoundaryMesh& bMesh = mesh_.boundaryMesh();

    labelHashSet visited;
    const label seedGId = globalNumbering_.toGlobal(seedCelli);

    // Queue of (localCellId, globalId). localCellId = -1 for remote.
    std::queue<std::pair<label, label>> bfsQueue;
    bfsQueue.push({seedCelli, seedGId});
    visited.insert(seedGId);

    // Seed at (0,0,0)
    topo[a3(0, 0, 0)] = seedGId;

    while (!bfsQueue.empty())
    {
        auto [localId, globalId] = bfsQueue.front();
        bfsQueue.pop();

        if (localId < 0) continue;  // remote: cannot expand

        const cell& c = mesh_.cells()[localId];

        forAll(c, fi)
        {
            const label faceI = c[fi];

            label nbrLocalId = -1;
            label nbrGlobalId = -1;
            point nbrCC;
            bool foundNbr = false;

            if (faceI < nInternalFaces)
            {
                nbrLocalId = (own[faceI] == localId)
                    ? nei[faceI] : own[faceI];
                nbrGlobalId = globalNumbering_.toGlobal(nbrLocalId);
                nbrCC = C_[nbrLocalId];
                foundNbr = true;
            }
            else
            {
                const label bFaceI = faceI - nInternalFaces;
                const label patchi = bMesh.patchID()[bFaceI];
                const polyPatch& pp = bMesh[patchi];

                if (isA<processorPolyPatch>(pp))
                {
                    const label localFaceI = faceI - pp.start();
                    if (procPatchRemoteGId_[patchi].size() > 0)
                    {
                        nbrGlobalId =
                            procPatchRemoteGId_[patchi][localFaceI];
                        nbrCC =
                            procPatchRemoteCC_[patchi][localFaceI];
                        nbrLocalId = -1;
                        foundNbr = true;
                    }
                }
                else if (
                    isA<cyclicPolyPatch>(pp)
                    && !isA<processorCyclicPolyPatch>(pp)
                )
                {
                    const cyclicPolyPatch& cyclicPatch =
                        refCast<const cyclicPolyPatch>(pp);
                    const cyclicPolyPatch& nbrPatch =
                        cyclicPatch.neighbPatch();
                    const label localFaceI = faceI - pp.start();
                    nbrLocalId = nbrPatch.faceCells()[localFaceI];
                    nbrGlobalId =
                        globalNumbering_.toGlobal(nbrLocalId);
                    nbrCC = C_[nbrLocalId];
                    foundNbr = true;
                }
                else if (isA<processorCyclicPolyPatch>(pp))
                {
                    const label localFaceI = faceI - pp.start();
                    if (procPatchRemoteGId_[patchi].size() > 0)
                    {
                        nbrGlobalId =
                            procPatchRemoteGId_[patchi][localFaceI];
                        nbrCC =
                            procPatchRemoteCC_[patchi][localFaceI];
                        nbrLocalId = -1;
                        foundNbr = true;
                    }
                }
            }

            if (!foundNbr) continue;
            if (visited.found(nbrGlobalId)) continue;
            visited.insert(nbrGlobalId);

            Vector<label> rel = relativeIJK(seedCC, nbrCC);
            if (!withinStencil(rel)) continue;

            const label flatIdx = a3(rel.x(), rel.y(), rel.z());
            topo[flatIdx] = nbrGlobalId;

            // Accumulate remote requests
            if (nbrLocalId < 0)
            {
                const label procID =
                    globalNumbering_.whichProcID(nbrGlobalId);
                commNeeded[procID].insert(nbrGlobalId);
            }

            // Continue BFS from local cells, stop at stencil corners
            if (nbrLocalId >= 0)
            {
                bool atCorner = (
                    (iMax_ > 0 && mag(rel.x()) == iMax_) &&
                    (jMax_ > 0 && mag(rel.y()) == jMax_) &&
                    (kMax_ > 0 && mag(rel.z()) == kMax_)
                );

                if (!atCorner)
                {
                    bfsQueue.push({nbrLocalId, nbrGlobalId});
                }
            }
            else
            {
                // Remote cell: add to BFS but mark as non-expandable
                bfsQueue.push({-1, nbrGlobalId});
            }
        }
    }

    return topo;
}


Foam::Vector<Foam::label> Foam::amrStencil::relativeIJK
(
    const point& seedCC,
    const point& targetCC
) const
{
    vector delta = targetCC - seedCC;

    if (isCyclic_.x() && Lx_ > SMALL)
    {
        if (delta.x() > 0.5 * Lx_)       delta.x() -= Lx_;
        else if (delta.x() < -0.5 * Lx_) delta.x() += Lx_;
    }
    if (isCyclic_.y() && Ly_ > SMALL)
    {
        if (delta.y() > 0.5 * Ly_)       delta.y() -= Ly_;
        else if (delta.y() < -0.5 * Ly_) delta.y() += Ly_;
    }
    if (isCyclic_.z() && Lz_ > SMALL)
    {
        if (delta.z() > 0.5 * Lz_)       delta.z() -= Lz_;
        else if (delta.z() < -0.5 * Lz_) delta.z() += Lz_;
    }

    label i = (isEmpty_.x()) ? 0 : round(delta.x() / dx_);
    label j = (isEmpty_.y()) ? 0 : round(delta.y() / dx_);
    label k = (isEmpty_.z()) ? 0 : round(delta.z() / dx_);

    return Vector<label>(i, j, k);
}


bool Foam::amrStencil::withinStencil(const Vector<label>& rel) const
{
    return (
        rel.x() >= -iMax_ && rel.x() <= iMax_ &&
        rel.y() >= -jMax_ && rel.y() <= jMax_ &&
        rel.z() >= -kMax_ && rel.z() <= kMax_
    );
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::amrStencil::updateIfNeeded()
{
    if (
        mesh_.nCells() != nCellsAtUpdate_
        || mesh_.time().timeIndex() != topoVersion_
            && mesh_.topoChanging()
    )
    {
        update();
    }
}


void Foam::amrStencil::update()
{
    // Clear all cached data
    stencilTopology_.clear();
    procMapsBuilt_ = false;

    // Rebuild global numbering for new mesh size
    globalNumbering_.reset(mesh_.nCells() + mesh_.nBoundaryFaces());

    // Rebuild neighbour proc list (topology may have changed)
    if (Pstream::parRun())
    {
        const polyBoundaryMesh& bMesh = mesh_.boundaryMesh();
        HashSet<label> procNeighbours;

        forAll(bMesh, patchi)
        {
            const polyPatch& pp = bMesh[patchi];
            if (isA<processorPolyPatch>(pp))
            {
                const processorPolyPatch& procPatch =
                    refCast<const processorPolyPatch>(pp);
                procNeighbours.insert(procPatch.neighbProcNo());
            }
        }
        neighbourProcs_ = procNeighbours.toc();
    }

    // Rebuild proc boundary face maps
    buildProcMaps();

    // Store current state
    topoVersion_ = mesh_.time().timeIndex();
    nCellsAtUpdate_ = mesh_.nCells();

    Info << "amrStencil::update(): caches rebuilt at time index "
         << topoVersion_ << ", nCells=" << nCellsAtUpdate_ << endl;
}


void Foam::amrStencil::setStencil
(
    const Map<scalar>& phiIJK,
    const label seedCelli
)
{
    A_ = 0.0;

    // Get or build topology for this seed cell
    const labelList& topo = stencilTopology_[seedCelli];

    // Fill stencil from phiIJK using cached topology
    forAll(topo, idx)
    {
        const label gId = topo[idx];
        if (gId >= 0 && phiIJK.found(gId))
        {
            A_[idx] = phiIJK[gId];
        }
    }
}


// ************************************************************************* //
