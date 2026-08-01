/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2026 Asim Onder
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
#include <queue>
#include <cmath>

Foam::amrStencil::amrStencil
(
    const fvMesh& mesh,
    int iMax, int jMax, int kMax,
    scalar dx, scalar dy, scalar dz,
    const Vector<bool>& isEmpty,
    const Vector<bool>& isCyclic,
    scalar Lx, scalar Ly, scalar Lz
)
:
    uniformStencil(SubclassTag{}, mesh, iMax, jMax, kMax, isEmpty, dx),
    mesh_(mesh),
    globalNumbering_(mesh.nCells() + mesh.nBoundaryFaces()),
    topoVersion_(-1),
    nCellsAtUpdate_(0),
    dy_(dy), dz_(dz),
    isCyclic_(isCyclic),
    Lx_(Lx), Ly_(Ly), Lz_(Lz),
    procMapsBuilt_(false),
    alphaHashInvBin_(0),
    hashNx_(0),
    hashNy_(0),
    rescueMissingStencils_(false)
{
    boundBox bb = mesh_.bounds();
    reduce(bb.min(), minOp<point>());
    reduce(bb.max(), maxOp<point>());
    if (isEmpty.y()) dy_ = bb.max().y() - bb.min().y();
    if (isEmpty.z()) dz_ = bb.max().z() - bb.min().z();
    domainMin_ = bb.min();

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

    buildProcMaps();
    nCellsAtUpdate_ = mesh_.nCells();
    topoVersion_ = mesh_.time().timeIndex();

    commAlgorithm_ = "allToAll";
    rescueMissingStencils_ = false;
    {
        IOobject dictHeader
        (
            "amrStencilDict",
            mesh_.time().system(),
            mesh_.time().db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        );
        if (dictHeader.typeHeaderOk<IOdictionary>(true))
        {
            IOdictionary dict(dictHeader);
            commAlgorithm_ =
                dict.lookupOrDefault<word>("commAlgorithm", "allToAll");
            rescueMissingStencils_ =
                dict.lookupOrDefault<bool>("rescueMissingStencils", false);
        }
    }
    Info<< "amrStencil: commAlgorithm = " << commAlgorithm_
        << ", rescueMissingStencils = "
        << (rescueMissingStencils_ ? "true" : "false") << endl;
    Info << "amrStencil: iMax=" << iMax_ << " jMax=" << jMax_
         << " kMax=" << kMax_ << " N=" << N_
         << " nNeighbourProcs=" << neighbourProcs_.size() << endl;
}


void Foam::amrStencil::buildProcMaps()
{
    if (!Pstream::parRun()) { procMapsBuilt_ = true; return; }

    const polyBoundaryMesh& bMesh = mesh_.boundaryMesh();
    procPatchRemoteGId_.setSize(bMesh.size());
    procPatchRemoteCC_.setSize(bMesh.size());

    PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);
    forAll(bMesh, patchi)
    {
        const polyPatch& pp = bMesh[patchi];
        if (!isA<processorPolyPatch>(pp)) continue;
        const processorPolyPatch& procPatch =
            refCast<const processorPolyPatch>(pp);
        const labelUList& faceCells = procPatch.faceCells();
        labelList gIds(faceCells.size());
        pointField ccs(faceCells.size());
        forAll(faceCells, fi)
        {
            gIds[fi] = globalNumbering_.toGlobal(faceCells[fi]);
            ccs[fi] = mesh_.C()[faceCells[fi]];
        }
        UOPstream toNbr(procPatch.neighbProcNo(), pBufs);
        toNbr << gIds << ccs;
    }
    pBufs.finishedSends();
    forAll(bMesh, patchi)
    {
        const polyPatch& pp = bMesh[patchi];
        if (!isA<processorPolyPatch>(pp))
        {
            procPatchRemoteGId_[patchi].clear();
            procPatchRemoteCC_[patchi].clear();
            continue;
        }
        const processorPolyPatch& procPatch =
            refCast<const processorPolyPatch>(pp);
        UIPstream fromNbr(procPatch.neighbProcNo(), pBufs);
        fromNbr >> procPatchRemoteGId_[patchi]
                >> procPatchRemoteCC_[patchi];
    }
    procMapsBuilt_ = true;
}


Foam::labelList Foam::amrStencil::buildStencilTopology
(
    const label seedCelli,
    List<labelHashSet>& commNeeded
)
{
    // BFS kept for reconstruction schemes only.
    labelList topo(N_, -1);
    const point seedCC = mesh_.C()[seedCelli];
    const labelUList& own = mesh_.faceOwner();
    const labelUList& nei = mesh_.faceNeighbour();
    const label nInternalFaces = mesh_.nInternalFaces();
    const polyBoundaryMesh& bMesh = mesh_.boundaryMesh();

    const scalar V = mesh_.V()[seedCelli];
    scalar localDx, localDy, localDz;
    if (isEmpty_.x() && !isEmpty_.y() && !isEmpty_.z())
    { scalar dl = Foam::sqrt(V/dx_); localDx=dx_; localDy=dl; localDz=dl; }
    else if (!isEmpty_.x() && isEmpty_.y() && !isEmpty_.z())
    { scalar dl = Foam::sqrt(V/dy_); localDx=dl; localDy=dy_; localDz=dl; }
    else if (!isEmpty_.x() && !isEmpty_.y() && isEmpty_.z())
    { scalar dl = Foam::sqrt(V/dz_); localDx=dl; localDy=dl; localDz=dz_; }
    else
    { scalar dl = Foam::cbrt(V); localDx=dl; localDy=dl; localDz=dl; }

    labelHashSet visited;
    const label seedGId = globalNumbering_.toGlobal(seedCelli);
    std::queue<label> bfsQueue;
    bfsQueue.push(seedCelli);
    visited.insert(seedGId);
    topo[a3(0,0,0)] = seedGId;

    while (!bfsQueue.empty())
    {
        const label localId = bfsQueue.front(); bfsQueue.pop();
        const cell& c = mesh_.cells()[localId];
        forAll(c, fi)
        {
            const label faceI = c[fi];
            label nbrLocal=-1, nbrGlobal=-1;
            point nbrCC(Zero);
            bool found=false;

            if (faceI < nInternalFaces)
            {
                nbrLocal = (own[faceI]==localId) ? nei[faceI] : own[faceI];
                nbrGlobal = globalNumbering_.toGlobal(nbrLocal);
                nbrCC = mesh_.C()[nbrLocal]; found=true;
            }
            else
            {
                const label bFaceI = faceI - nInternalFaces;
                const label patchi = bMesh.patchID()[bFaceI];
                const polyPatch& pp = bMesh[patchi];
                if (isA<processorPolyPatch>(pp))
                {
                    const label lfi = faceI - pp.start();
                    if (procPatchRemoteGId_[patchi].size()>0)
                    { nbrGlobal=procPatchRemoteGId_[patchi][lfi];
                      nbrCC=procPatchRemoteCC_[patchi][lfi];
                      nbrLocal=-1; found=true; }
                }
                else if (isA<cyclicPolyPatch>(pp) && !isA<processorCyclicPolyPatch>(pp))
                {
                    const cyclicPolyPatch& cyc = refCast<const cyclicPolyPatch>(pp);
                    const label lfi = faceI - pp.start();
                    nbrLocal = cyc.neighbPatch().faceCells()[lfi];
                    nbrGlobal = globalNumbering_.toGlobal(nbrLocal);
                    nbrCC = mesh_.C()[nbrLocal]; found=true;
                }
                else if (isA<processorCyclicPolyPatch>(pp))
                {
                    const label lfi = faceI - pp.start();
                    if (procPatchRemoteGId_[patchi].size()>0)
                    { nbrGlobal=procPatchRemoteGId_[patchi][lfi];
                      nbrCC=procPatchRemoteCC_[patchi][lfi];
                      nbrLocal=-1; found=true; }
                }
            }
            if (!found) continue;
            if (visited.found(nbrGlobal)) continue;
            visited.insert(nbrGlobal);

            vector delta = nbrCC - seedCC;
            if (isCyclic_.x()&&Lx_>SMALL)
            { if(delta.x()>0.5*Lx_) delta.x()-=Lx_;
              else if(delta.x()<-0.5*Lx_) delta.x()+=Lx_; }
            if (isCyclic_.y()&&Ly_>SMALL)
            { if(delta.y()>0.5*Ly_) delta.y()-=Ly_;
              else if(delta.y()<-0.5*Ly_) delta.y()+=Ly_; }
            if (isCyclic_.z()&&Lz_>SMALL)
            { if(delta.z()>0.5*Lz_) delta.z()-=Lz_;
              else if(delta.z()<-0.5*Lz_) delta.z()+=Lz_; }

            Vector<label> rel
            (
                isEmpty_.x()?0:round(delta.x()/localDx),
                isEmpty_.y()?0:round(delta.y()/localDy),
                isEmpty_.z()?0:round(delta.z()/localDz)
            );
            if (!withinStencil(rel)) continue;
            topo[a3(rel.x(),rel.y(),rel.z())] = nbrGlobal;
            if (!globalNumbering_.isLocal(nbrGlobal))
                commNeeded[globalNumbering_.whichProcID(nbrGlobal)].insert(nbrGlobal);
            if (nbrLocal>=0) bfsQueue.push(nbrLocal);
        }
    }
    return topo;
}


bool Foam::amrStencil::withinStencil(const Vector<label>& rel) const
{
    return (rel.x()>=-iMax_ && rel.x()<=iMax_ &&
            rel.y()>=-jMax_ && rel.y()<=jMax_ &&
            rel.z()>=-kMax_ && rel.z()<=kMax_);
}


void Foam::amrStencil::updateIfNeeded()
{
    if (mesh_.nCells() != nCellsAtUpdate_
        || (mesh_.time().timeIndex() != topoVersion_ && mesh_.topoChanging()))
    {
        update();
    }
}


void Foam::amrStencil::update()
{
    stencilTopology_.clear();
    sendCellsPerProc_.clear();
    exchangeProcs_.clear();
    interfaceBandCells_.clear();
    procMapsBuilt_ = false;

    // Clear alpha hash
    alphaHash_.clear();
    alphaHashInvBin_ = 0;
    hashNx_ = 0;
    hashNy_ = 0;

    globalNumbering_.reset(mesh_.nCells() + mesh_.nBoundaryFaces());

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
    buildProcMaps();
    topoVersion_ = mesh_.time().timeIndex();
    nCellsAtUpdate_ = mesh_.nCells();
}


Foam::point Foam::amrStencil::cyclicWrap(const point& p) const
{
    point q = p;
    const point& bMin = domainMin_;
    if (isCyclic_.x() && Lx_ > SMALL)
    {
        if (q.x() < bMin.x()) q.x() += Lx_;
        else if (q.x() > bMin.x() + Lx_) q.x() -= Lx_;
    }
    if (isCyclic_.y() && Ly_ > SMALL)
    {
        if (q.y() < bMin.y()) q.y() += Ly_;
        else if (q.y() > bMin.y() + Ly_) q.y() -= Ly_;
    }
    if (isCyclic_.z() && Lz_ > SMALL)
    {
        if (q.z() < bMin.z()) q.z() += Lz_;
        else if (q.z() > bMin.z() + Lz_) q.z() -= Lz_;
    }
    return q;
}


void Foam::amrStencil::setStencil
(
    const Map<scalar>& phiIJK,
    const label seedCelli
)
{
    if (!stencilTopology_.found(seedCelli))
    {
        FatalErrorInFunction
            << "No topology for cell " << seedCelli
            << abort(FatalError);
    }
    const labelList& topo = stencilTopology_[seedCelli];
    forAll(topo, idx)
    {
        const label gId = topo[idx];
        if (gId >= 0 && phiIJK.found(gId))
            A_[idx] = phiIJK[gId];
        else
            A_[idx] = 0;
    }
}

// ************************************************************************* //
