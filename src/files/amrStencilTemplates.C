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
    // Auto-detect mesh changes and rebuild caches if needed
    updateIfNeeded();

    // Ensure proc maps are built
    if (!procMapsBuilt_)
    {
        buildProcMaps();
    }

    List<labelHashSet> commNeeded(Pstream::nProcs());

    // ---------------------------------------------------------------
    // Phase 1: For each interface cell, get or build stencil topology.
    //          Collect local field values directly, accumulate remote
    //          cell requests.
    // ---------------------------------------------------------------

    forAll(zone, celli)
    {
        if (!zone[celli]) continue;

        // Build and cache topology if not already cached
        if (!stencilTopology_.found(celli))
        {
            stencilTopology_.insert
            (
                celli,
                buildStencilTopology(celli, commNeeded)
            );
        }
        else
        {
            // Topology is cached. Still need to check for remote cells
            // that need to be requested. Scan the cached topology.
            const labelList& topo = stencilTopology_[celli];

            forAll(topo, idx)
            {
                const label gId = topo[idx];
                if (gId < 0) continue;

                if (!globalNumbering_.isLocal(gId))
                {
                    const label procID =
                        globalNumbering_.whichProcID(gId);
                    commNeeded[procID].insert(gId);
                }
            }
        }

        // Collect local values from the topology
        const labelList& topo = stencilTopology_[celli];

        forAll(topo, idx)
        {
            const label gId = topo[idx];
            if (gId < 0) continue;

            if (globalNumbering_.isLocal(gId) && !phiMap.found(gId))
            {
                phiMap.insert
                (
                    gId,
                    phi[globalNumbering_.toLocal(gId)]
                );
            }
        }
    }

    // ---------------------------------------------------------------
    // Phase 2: Parallel exchange (same pattern as zoneTransfers.H)
    // ---------------------------------------------------------------

    if (Pstream::parRun())
    {
        List<labelHashSet> send(Pstream::nProcs());

        // Tell each neighbour proc which global IDs we need
        {
            PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

            forAll(neighbourProcs_, i)
            {
                const label domain = neighbourProcs_[i];
                UOPstream toDomain(domain, pBufs);
                toDomain << commNeeded[domain];
            }

            pBufs.finishedSends();

            forAll(neighbourProcs_, i)
            {
                const label domain = neighbourProcs_[i];
                send[domain].clear();
                UIPstream fromDomain(domain, pBufs);
                fromDomain >> send[domain];
            }
        }

        // Prepare and send field values
        List<Map<Type>> sendValues(Pstream::nProcs());

        forAll(neighbourProcs_, i)
        {
            const label domaini = neighbourProcs_[i];
            for (const label sendIdx : send[domaini])
            {
                sendValues[domaini].insert
                (
                    sendIdx,
                    phi[globalNumbering_.toLocal(sendIdx)]
                );
            }
        }

        {
            PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

            forAll(neighbourProcs_, i)
            {
                const label domain = neighbourProcs_[i];
                UOPstream toDomain(domain, pBufs);
                toDomain << sendValues[domain];
            }

            pBufs.finishedSends();

            forAll(neighbourProcs_, i)
            {
                const label domain = neighbourProcs_[i];
                Map<Type> tmpValue;
                UIPstream fromDomain(domain, pBufs);
                fromDomain >> tmpValue;
                phiMap += tmpValue;
            }
        }
    }
}


// ************************************************************************* //
