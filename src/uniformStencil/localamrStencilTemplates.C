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

#include <queue>

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<typename Type>
void Foam::amrStencil::getZoneField
(
    const boolList& zone,
    Map<Type>& phiMap,
    const GeometricField<Type, fvPatchField, volMesh>& phi
)
{
    updateIfNeeded();

    if (!procMapsBuilt_)
    {
        buildProcMaps();
    }

    const label nLayers = max(max(iMax_, jMax_), kMax_);

    // ---------------------------------------------------------------
    // Step 0 (parallel): Scalable exchange.
    //   (a) Gather all procs' interface bboxes (2 reduces).
    //   (b) Each proc checks which remote interface bboxes are
    //       within stencil reach of its OWN domain -- those procs
    //       might need my cells.  Symmetrise via allGather.
    //   (c) For each partner, send only cells within that partner's
    //       interface bbox expanded by stencil reach.
    // ---------------------------------------------------------------

    labelList  remoteGIds;
    pointField remoteCCs;
    Map<Type>  gIdToValue;

    if (Pstream::parRun())
    {
        // --- Gather interface bounding boxes from all procs ---

        point localIMin = point::uniform(VGREAT);
        point localIMax = point::uniform(-VGREAT);

        forAll(zone, ci)
        {
            if (zone[ci])
            {
                const point& cc = mesh_.C()[ci];
                localIMin = Foam::min(localIMin, cc);
                localIMax = Foam::max(localIMax, cc);
            }
        }

        List<point> allIMin(Pstream::nProcs(), point::uniform(VGREAT));
        List<point> allIMax(Pstream::nProcs(), point::uniform(-VGREAT));

        allIMin[Pstream::myProcNo()] = localIMin;
        allIMax[Pstream::myProcNo()] = localIMax;

        Foam::reduce(allIMin, minOp<List<point>>());
        Foam::reduce(allIMax, maxOp<List<point>>());

        // --- Estimate cell spacing ---

        scalar dl = 0;
        forAll(zone, ci)
        {
            if (zone[ci])
            {
                const scalar V = mesh_.V()[ci];
                dl = (isEmpty_.x() || isEmpty_.y() || isEmpty_.z())
                    ? Foam::sqrt(V / max(dx_, max(dy_, dz_)))
                    : Foam::cbrt(V);
                break;
            }
        }
        reduce(dl, maxOp<scalar>());

        // Stencil reach: nLayers cells in each direction.
        // Use 1.5x safety margin to account for AMR cell size variation.
        const scalar reach = 1.5 * nLayers * dl;
        const vector L(Lx_, Ly_, Lz_);

        // --- My domain bounding box ---

        const boundBox myBB = mesh_.bounds();

        // --- Determine exchange partners ---
        // Proc B is my partner if:
        //   - B has interface cells, AND
        //   - B's interface bbox expanded by `reach` overlaps my domain
        //     (B's stencils might need my cells), OR
        //   - I have interface cells AND my interface bbox expanded by
        //     `reach` overlaps B's domain (I might need B's cells).
        //
        // We compute one side (who might need MY cells) and then
        // symmetrise via boolList allGather so both sides agree.

        // Step 1: Which remote procs' interfaces are within reach
        // of my domain?  Those procs might need my cells.
        boolList iWantToSendTo(Pstream::nProcs(), false);

        for (label p = 0; p < Pstream::nProcs(); p++)
        {
            if (p == Pstream::myProcNo()) continue;
            if (allIMin[p].x() > VGREAT/2) continue;  // no interface

            // Expand remote interface bbox by reach
            point eMin = allIMin[p] - point::uniform(reach);
            point eMax = allIMax[p] + point::uniform(reach);

            // Check overlap with my domain bbox, accounting for cyclics
            bool overlaps = true;
            for (direction d = 0; d < 3; d++)
            {
                bool dimOverlap = false;

                // Direct overlap
                if (eMin[d] <= myBB.max()[d]
                    && eMax[d] >= myBB.min()[d])
                {
                    dimOverlap = true;
                }

                // Cyclic shifted overlap
                if (!dimOverlap && isCyclic_[d] && L[d] > SMALL)
                {
                    // Shift expanded bbox by +L
                    if ((eMin[d] + L[d]) <= myBB.max()[d]
                        && (eMax[d] + L[d]) >= myBB.min()[d])
                    {
                        dimOverlap = true;
                    }
                    // Shift expanded bbox by -L
                    if ((eMin[d] - L[d]) <= myBB.max()[d]
                        && (eMax[d] - L[d]) >= myBB.min()[d])
                    {
                        dimOverlap = true;
                    }
                }

                if (!dimOverlap)
                {
                    overlaps = false;
                    break;
                }
            }

            if (overlaps) iWantToSendTo[p] = true;
        }

        // Step 2: Symmetrise.  Gather all procs' send intentions.
        // wantsMatrix[p] = boolList of size nProcs, true if p wants
        // to send to that proc.
        // Then: I exchange with proc p if I want to send to p OR
        // p wants to send to me.

        List<boolList> wantsMatrix(Pstream::nProcs());
        wantsMatrix[Pstream::myProcNo()] = iWantToSendTo;

        // AllGather via PstreamBuffers
        {
            PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

            for (label p = 0; p < Pstream::nProcs(); p++)
            {
                if (p != Pstream::myProcNo())
                {
                    UOPstream os(p, pBufs);
                    os << iWantToSendTo;
                }
            }

            pBufs.finishedSends();

            for (label p = 0; p < Pstream::nProcs(); p++)
            {
                if (p != Pstream::myProcNo())
                {
                    UIPstream is(p, pBufs);
                    is >> wantsMatrix[p];
                }
            }
        }

        HashSet<label> exchangeSet;
        for (label p = 0; p < Pstream::nProcs(); p++)
        {
            if (p == Pstream::myProcNo()) continue;
            // I exchange with p if I want to send to p,
            // OR p wants to send to me.
            if (iWantToSendTo[p]
                || wantsMatrix[p][Pstream::myProcNo()])
            {
                exchangeSet.insert(p);
            }
        }

        exchangeProcs_ = exchangeSet.toc();

        Pout << "amrStencil: proc " << Pstream::myProcNo()
             << " nLocalCells=" << mesh_.nCells()
             << " nExchangeProcs=" << exchangeProcs_.size()
             << " dl=" << dl << " reach=" << reach
             << endl;

        // ----------------------------------------------------------
        // Step 0b: Exchange cells.
        //   For each partner p, send local cells that fall within
        //   p's interface bbox expanded by `reach`.
        //   If p has no interface (it's a pure sender), send nothing.
        // ----------------------------------------------------------

        {
            PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

            const pointField& cc = mesh_.C();

            forAll(exchangeProcs_, ni)
            {
                const label p = exchangeProcs_[ni];

                // Determine which local cells to send to proc p.
                // If proc p has interface cells, send cells within
                // its expanded interface bbox.
                // If proc p has no interface, it won't need our cells
                // (it only sends), so send empty.

                labelList sGIds;
                pointField sCCs;
                List<Type> sVals;

                if (allIMin[p].x() < VGREAT/2)  // p has interface
                {
                    point eMin = allIMin[p] - point::uniform(reach);
                    point eMax = allIMax[p] + point::uniform(reach);

                    DynamicList<label> sendGIds;
                    DynamicList<point> sendCCs;
                    DynamicList<Type>  sendVals;

                    for (label ci = 0; ci < mesh_.nCells(); ci++)
                    {
                        bool inside = true;
                        for (direction d = 0; d < 3; d++)
                        {
                            bool dimInside = false;

                            // Direct
                            if (cc[ci][d] >= eMin[d]
                                && cc[ci][d] <= eMax[d])
                            {
                                dimInside = true;
                            }

                            // Cyclic shifted
                            if (!dimInside && isCyclic_[d]
                                && L[d] > SMALL)
                            {
                                scalar shifted = cc[ci][d] + L[d];
                                if (shifted >= eMin[d]
                                    && shifted <= eMax[d])
                                {
                                    dimInside = true;
                                }
                                shifted = cc[ci][d] - L[d];
                                if (shifted >= eMin[d]
                                    && shifted <= eMax[d])
                                {
                                    dimInside = true;
                                }
                            }

                            if (!dimInside)
                            {
                                inside = false;
                                break;
                            }
                        }

                        if (inside)
                        {
                            sendGIds.append(
                                globalNumbering_.toGlobal(ci));
                            sendCCs.append(cc[ci]);
                            sendVals.append(phi[ci]);
                        }
                    }

                    sGIds = sendGIds;
                    sCCs = sendCCs;
                    sVals = sendVals;
                }

                UOPstream toNbr(p, pBufs);
                toNbr << sGIds << sCCs << sVals;
            }

            pBufs.finishedSends();

            // Receive
            DynamicList<label> allGIds;
            DynamicList<point> allCCs;

            forAll(exchangeProcs_, ni)
            {
                labelList recvG;
                pointField recvC;
                List<Type> recvV;

                UIPstream fromNbr(exchangeProcs_[ni], pBufs);
                fromNbr >> recvG >> recvC >> recvV;

                forAll(recvG, ci)
                {
                    const label gId = recvG[ci];
                    if (!phiMap.found(gId))
                    {
                        phiMap.insert(gId, recvV[ci]);
                    }
                    if (!gIdToValue.found(gId))
                    {
                        gIdToValue.insert(gId, recvV[ci]);
                        allGIds.append(gId);
                        allCCs.append(recvC[ci]);
                    }
                }
            }

            remoteGIds = allGIds;
            remoteCCs = allCCs;
        }
    }

    // ---------------------------------------------------------------
    // Step 1: Build topology (local BFS), collect local values,
    //         nearest-neighbour fill for remote cells.
    // ---------------------------------------------------------------

    List<labelHashSet> commDummy(Pstream::nProcs());

    forAll(zone, celli)
    {
        if (!zone[celli]) continue;

        if (!stencilTopology_.found(celli))
        {
            stencilTopology_.insert
            (
                celli,
                buildStencilTopology(celli, commDummy)
            );
        }

        const labelList& topo = stencilTopology_[celli];

        // Collect local values
        forAll(topo, idx)
        {
            const label gId = topo[idx];
            if (gId >= 0 && globalNumbering_.isLocal(gId)
                && !phiMap.found(gId))
            {
                phiMap.insert(gId, phi[globalNumbering_.toLocal(gId)]);
            }
        }

        // Nearest-neighbour fill for remote/missing stencil positions
        if (Pstream::parRun() && remoteCCs.size() > 0)
        {
            const point seedCC = mesh_.C()[celli];
            const scalar V = mesh_.V()[celli];
            scalar localDx, localDy, localDz;

            if (isEmpty_.x() && !isEmpty_.y() && !isEmpty_.z())
            {
                scalar dl = Foam::sqrt(V / dx_);
                localDx = dx_; localDy = dl; localDz = dl;
            }
            else if (!isEmpty_.x() && isEmpty_.y() && !isEmpty_.z())
            {
                scalar dl = Foam::sqrt(V / dy_);
                localDx = dl; localDy = dy_; localDz = dl;
            }
            else if (!isEmpty_.x() && !isEmpty_.y() && isEmpty_.z())
            {
                scalar dl = Foam::sqrt(V / dz_);
                localDx = dl; localDy = dl; localDz = dz_;
            }
            else
            {
                scalar dl = Foam::cbrt(V);
                localDx = dl; localDy = dl; localDz = dl;
            }

            const scalar tol = 1.5 * min(min(localDx, localDy), localDz);
            const scalar tolSq = tol * tol;

            for (label i = -iMax_; i <= iMax_; i++)
            for (label j = -jMax_; j <= jMax_; j++)
            for (label k = -kMax_; k <= kMax_; k++)
            {
                const label flatIdx = a3(i, j, k);
                const label gId = topo[flatIdx];

                bool needFill = (gId < 0);
                if (!needFill && !globalNumbering_.isLocal(gId)
                    && !phiMap.found(gId))
                {
                    needFill = true;
                }
                if (!needFill) continue;

                point expectedCC = seedCC;
                if (!isEmpty_.x()) expectedCC.x() += i * localDx;
                if (!isEmpty_.y()) expectedCC.y() += j * localDy;
                if (!isEmpty_.z()) expectedCC.z() += k * localDz;

                // Cyclic wrapping
                const point& bMin = domainMin_;
                if (isCyclic_.x() && Lx_ > SMALL)
                {
                    if (expectedCC.x() < bMin.x())
                        expectedCC.x() += Lx_;
                    else if (expectedCC.x() > bMin.x() + Lx_)
                        expectedCC.x() -= Lx_;
                }
                if (isCyclic_.y() && Ly_ > SMALL)
                {
                    if (expectedCC.y() < bMin.y())
                        expectedCC.y() += Ly_;
                    else if (expectedCC.y() > bMin.y() + Ly_)
                        expectedCC.y() -= Ly_;
                }
                if (isCyclic_.z() && Lz_ > SMALL)
                {
                    if (expectedCC.z() < bMin.z())
                        expectedCC.z() += Lz_;
                    else if (expectedCC.z() > bMin.z() + Lz_)
                        expectedCC.z() -= Lz_;
                }

                // Nearest-neighbour search
                scalar bestDistSq = GREAT;
                label  bestIdx = -1;

                forAll(remoteCCs, ri)
                {
                    const scalar dSq =
                        Foam::magSqr(remoteCCs[ri] - expectedCC);
                    if (dSq < bestDistSq)
                    {
                        bestDistSq = dSq;
                        bestIdx = ri;
                    }
                }

                if (bestIdx >= 0 && bestDistSq < tolSq)
                {
                    const label foundGId = remoteGIds[bestIdx];
                    if (!phiMap.found(foundGId)
                        && gIdToValue.found(foundGId))
                    {
                        phiMap.insert(foundGId, gIdToValue[foundGId]);
                    }
                    stencilTopology_[celli][flatIdx] = foundGId;
                }
            }
        }
    }
}


// ************************************************************************* //
