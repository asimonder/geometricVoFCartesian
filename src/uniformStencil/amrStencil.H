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

Class
    Foam::amrStencil

Description
    Stencil class for height function curvature on AMR meshes.
    Builds a direct position-to-alpha hash for O(1) stencil lookups.
    Remote cell values obtained via face-adjacent exchange with
    round-2 relay.

SourceFiles
    amrStencil.C
    amrStencilTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef amrStencil_H
#define amrStencil_H

#include "uniformStencil.H"
#include "globalIndex.H"
#include "processorPolyPatch.H"
#include "processorCyclicPolyPatch.H"
#include "cyclicPolyPatch.H"
#include <cmath>

namespace Foam
{

class amrStencil
:
    public uniformStencil
{
private:

    const fvMesh& mesh_;

    globalIndex globalNumbering_;

    label topoVersion_;
    label nCellsAtUpdate_;

    scalar dy_;
    scalar dz_;

    Vector<bool> isCyclic_;
    scalar Lx_, Ly_, Lz_;
    point domainMin_;

    labelList neighbourProcs_;

    // Processor boundary data
    List<labelList> procPatchRemoteGId_;
    List<pointField> procPatchRemoteCC_;
    bool procMapsBuilt_;

    // Direct position-to-alpha hash
    // Key: floor(cc.x*invBin) + 99991*floor(cc.y*invBin) + 999979*floor(cc.z*invBin)
    // Value: alpha at that cell
    Map<scalar> alphaHash_;
    scalar alphaHashInvBin_;
    label hashNx_;
    label hashNy_;

    // Legacy stencil topology (kept for reconstruction schemes)
    Map<labelList> stencilTopology_;

    // Cached exchange data
    Map<labelList> sendCellsPerProc_;
    labelList exchangeProcs_;
    labelList interfaceBandCells_;

    word commAlgorithm_;
    bool rescueMissingStencils_;

    // Private Member Functions
    void buildProcMaps();
    void updateIfNeeded();

    labelList buildStencilTopology
    (
        const label seedCelli,
        List<labelHashSet>& commNeeded
    );

    bool withinStencil(const Vector<label>& rel) const;

    amrStencil(const amrStencil&) = delete;
    void operator=(const amrStencil&) = delete;

public:

    amrStencil
    (
        const fvMesh& mesh,
        int iMax, int jMax, int kMax,
        scalar dx, scalar dy, scalar dz,
        const Vector<bool>& isEmpty,
        const Vector<bool>& isCyclic,
        scalar Lx, scalar Ly, scalar Lz
    );

    void update();

    //- True when a stencil topology is cached for this seed cell
    bool hasStencil(const label seedCelli) const
    {
        return stencilTopology_.found(seedCelli);
    }

    void setStencil(const Map<scalar>& phiMap, const label seedCelli);

    template<typename Type>
    void getZoneField
    (
        const boolList& zone,
        Map<Type>& phiMap,
        const GeometricField<Type, fvPatchField, volMesh>& phi
    );

    //- Apply cyclic wrapping to a point
    point cyclicWrap(const point& p) const;

    //- Access the alpha hash for direct lookup
    const Map<scalar>& alphaHash() const { return alphaHash_; }

    //- Access inverse bin size for key computation
    scalar alphaHashInvBin() const { return alphaHashInvBin_; }

    //- Fill the stencil cube from the spatial hash around a seed cell
    //  centre with spacing h; empty directions stay at offset zero.
    //  Missing entries are set to zero; returns their count so the
    //  caller can fall back to another normal.
    label setStencilFromHash(const point& seedCC, const scalar h)
    {
        label nMissing = 0;
        for (label k = -kMax_; k <= kMax_; ++k)
        {
            for (label j = -jMax_; j <= jMax_; ++j)
            {
                for (label i = -iMax_; i <= iMax_; ++i)
                {
                    point p = seedCC + vector(i*h, j*h, k*h);
                    p = cyclicWrap(p);
                    const label key = posKey(p);
                    const label idx = a3(i, j, k);
                    const auto iter = alphaHash_.cfind(key);
                    if (iter.good())
                    {
                        A_[idx] = iter.val();
                    }
                    else
                    {
                        A_[idx] = 0;
                        ++nMissing;
                    }
                }
            }
        }
        return nMissing;
    }

    //- Compute hash key for a position (structured, collision-free)
    //  key = ix + Nx*iy + Nx*Ny*iz. Sequential, small keys.
    label posKey
    (
        const point& p
    ) const
    {
        const label ix = label(std::floor(
            (p.x() - domainMin_.x()) * alphaHashInvBin_));
        const label iy = label(std::floor(
            (p.y() - domainMin_.y()) * alphaHashInvBin_));
        const label iz = label(std::floor(
            (p.z() - domainMin_.z()) * alphaHashInvBin_));
        return ix + hashNx_ * (iy + hashNy_ * iz);
    }

    //- Access domain origin
    const point& domainMin() const { return domainMin_; }

    const labelList& getTopology(const label celli) const
    {
        return stencilTopology_[celli];
    }

    const globalIndex& globalNumbering() const
    {
        return globalNumbering_;
    }

    label globalNum() const
    {
        return globalNumbering_.totalSize();
    }
};

} // End namespace Foam

#ifdef NoRepository
    #include "amrStencilTemplates.C"
#endif

#endif
