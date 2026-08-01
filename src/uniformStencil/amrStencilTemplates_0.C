/*---------------------------------------------------------------------------* \
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

    // ===============================================================
    //  PARALLEL EXCHANGE
    //
    //  Partner list: ALL other procs (same as brute force that works).
    //  Send list: filtered per partner by partner's interface bbox
    //  expanded by stencil reach.
    //
    //  If this version works: partner detection was the bug.
    //  If this version fails: cell selection (bbox filter) is the bug.
    // ===============================================================

    labelList  remoteGIds;
    pointField remoteCCs;
    Map<Type>  gIdToValue;

    if (Pstream::parRun())
    {
        // --- Gather interface bboxes from all procs ---

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

        // Stencil reach: nLayers cells
        const scalar reach = nLayers * dl;
        const vector L(Lx_, Ly_, Lz_);

        // --- Build GLOBAL interface bbox (union of all procs) ---
        // A stencil can extend `reach` from any interface cell.
        // So any cell within `reach` of the GLOBAL interface bbox
        // could be needed by some proc.

        point globalIMin = point::uniform(VGREAT);
        point globalIMax = point::uniform(-VGREAT);
        for (label p = 0; p < Pstream::nProcs(); p++)
        {
            if (allIMin[p].x() < VGREAT/2)
            {
                globalIMin = Foam::min(globalIMin, allIMin[p]);
                globalIMax = Foam::max(globalIMax, allIMax[p]);
            }
        }

        // Expand by stencil reach
        point sendMin = globalIMin - point::uniform(reach);
        point sendMax = globalIMax + point::uniform(reach);

        // --- Partner list: ALL other procs ---
        labelList allOtherProcs(Pstream::nProcs() - 1);
        {
            label idx = 0;
            for (label p = 0; p < Pstream::nProcs(); p++)
            {
                if (p != Pstream::myProcNo())
                    allOtherProcs[idx++] = p;
            }
        }

        // --- Filter local cells: only those within expanded
        //     global interface bbox ---
        // This is the SAME set sent to every partner.
        // (Per-partner filtering is the next optimization step.)

        const pointField& cc = mesh_.C();

        DynamicList<label> sendGIds;
        DynamicList<point> sendCCs;
        DynamicList<Type>  sendVals;

        for (label ci = 0; ci < mesh_.nCells(); ci++)
        {
            bool inside = true;
            for (direction d = 0; d < 3; d++)
            {
                bool dimIn = false;

                // Direct
                if (cc[ci][d] >= sendMin[d]
                    && cc[ci][d] <= sendMax[d])
                {
                    dimIn = true;
                }

                // Cyclic shifted
                if (!dimIn && isCyclic_[d] && L[d] > SMALL)
                {
                    if ((cc[ci][d]+L[d]) >= sendMin[d]
                        && (cc[ci][d]+L[d]) <= sendMax[d])
                        dimIn = true;
                    if ((cc[ci][d]-L[d]) >= sendMin[d]
                        && (cc[ci][d]-L[d]) <= sendMax[d])
                        dimIn = true;
                }

                if (!dimIn)
                {
                    inside = false;
                    break;
                }
            }

            if (inside)
            {
                sendGIds.append(globalNumbering_.toGlobal(ci));
                sendCCs.append(cc[ci]);
                sendVals.append(phi[ci]);
            }
        }

        Pout << "amrStencil: proc " << Pstream::myProcNo()
             << " nLocalCells=" << mesh_.nCells()
             << " nSendCells=" << sendGIds.size()
             << " dl=" << dl << " reach=" << reach
             << " sendBox=(" << sendMin << " " << sendMax << ")"
             << endl;

        // --- Exchange ---
        {
            PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

            labelList sG(sendGIds);
            pointField sC(sendCCs);
            List<Type> sV(sendVals);

            forAll(allOtherProcs, ni)
            {
                UOPstream toNbr(allOtherProcs[ni], pBufs);
                toNbr << sG << sC << sV;
            }

            pBufs.finishedSends();

            DynamicList<label> allGIds;
            DynamicList<point> allCCs;

            forAll(allOtherProcs, ni)
            {
                labelList recvG;
                pointField recvC;
                List<Type> recvV;

                UIPstream fromNbr(allOtherProcs[ni], pBufs);
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
    // Step 1: Build topology, collect local values,
    //         nearest-neighbour fill for remote cells.
    // ---------------------------------------------------------------

    List<labelHashSet> commDummy(Pstream::nProcs());

    label nIncomplete = 0;

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

            label nMissing = 0;

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
                else
                {
                    nMissing++;
                    // Print first few misses with detail
                    if (nMissing <= 3 && nIncomplete < 5)
                    {
                        Pout << "  MISS cell=" << celli
                             << " ijk=(" << i << " " << j << " " << k
                             << ") expected=" << expectedCC
                             << " bestDist="
                             << (bestIdx >= 0 ? Foam::sqrt(bestDistSq) : -1)
                             << " tol=" << tol
                             << " gId=" << gId
                             << endl;
                    }
                }
            }

            if (nMissing > 0)
            {
                nIncomplete++;
            }
        }
    }

    if (nIncomplete > 0)
    {
        Pout << "amrStencil: proc " << Pstream::myProcNo()
             << " nIncompleteStencils=" << nIncomplete
             << endl;
    }
}


// ************************************************************************* //
