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

    // ---------------------------------------------------------------
    // Step 0 (parallel): Every proc sends ALL its cells to every
    //         other proc.  Brute force, guaranteed correct.
    //         Optimize later after correctness is verified.
    // ---------------------------------------------------------------

    labelList  remoteGIds;
    pointField remoteCCs;
    Map<Type>  gIdToValue;

    if (Pstream::parRun())
    {
        // Build list of ALL other procs
        labelList allOtherProcs(Pstream::nProcs() - 1);
        {
            label idx = 0;
            for (label p = 0; p < Pstream::nProcs(); p++)
            {
                if (p != Pstream::myProcNo())
                {
                    allOtherProcs[idx++] = p;
                }
            }
        }

        Pout << "amrStencil: proc " << Pstream::myProcNo()
             << " nLocalCells=" << mesh_.nCells()
             << " exchanging with ALL " << allOtherProcs.size()
             << " procs" << endl;

        // Send ALL local cells to every other proc
        PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

        const label nLocal = mesh_.nCells();
        labelList sGIds(nLocal);
        pointField sCCs(nLocal);
        List<Type> sVals(nLocal);

        for (label ci = 0; ci < nLocal; ci++)
        {
            sGIds[ci] = globalNumbering_.toGlobal(ci);
            sCCs[ci] = mesh_.C()[ci];
            sVals[ci] = phi[ci];
        }

        forAll(allOtherProcs, ni)
        {
            UOPstream toNbr(allOtherProcs[ni], pBufs);
            toNbr << sGIds << sCCs << sVals;
        }

        pBufs.finishedSends();

        // Receive from every other proc
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
