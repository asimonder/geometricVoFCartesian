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

Class
    Foam::amrStencil

Description
    AMR-compatible extension of uniformStencil. Discovers stencil cells
    via BFS neighbor traversal instead of global (i,j,k) lookup. Assigns
    relative (i,j,k) offsets geometrically, which is valid because the
    interfacial band is guaranteed uniform under AMR.

    Inherits all stencil operations (A_, a2, a3, calcYoungNormal,
    calcCCDNormal, flip, transpose, getSmallStencil) from uniformStencil.

    Caching strategy:
      - Processor boundary maps and per-cell stencil topology are cached.
      - Call update() after mesh changes (AMR refine/coarsen) to rebuild.
      - Between mesh changes, setStencil and getZoneField reuse cached
        topology with zero BFS cost.

    Developed by Asim Onder (2025)

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
#include "PstreamBuffers.H"
#include "HashSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class amrStencil Declaration
\*---------------------------------------------------------------------------*/

class amrStencil
:
    public uniformStencil
{
private:

    // Private data

    //- Reference to mesh
    const fvMesh& mesh_;

    //- Cell center field
    const volVectorField& C_;

    //- Global cell numbering for parallel lookup (rebuilt on mesh change)
    globalIndex globalNumbering_;

    //- Mesh topology version at last update. Used to auto-detect mesh
    //  changes and trigger cache invalidation.
    label topoVersion_;

    //- Number of cells at last update (for AMR detection)
    label nCellsAtUpdate_;

    //- Cyclic direction flags
    Vector<bool> isCyclic_;

    //- Domain extents (for cyclic wrapping)
    scalar Lx_;
    scalar Ly_;
    scalar Lz_;

    //- List of neighbour processor IDs
    labelList neighbourProcs_;


    // -- Cached processor boundary data (rebuilt on update()) --

    //- Face-indexed remote global IDs per processor patch
    //  procPatchRemoteGId_[patchi][localFaceI] = remote global ID
    List<labelList> procPatchRemoteGId_;

    //- Face-indexed remote cell centres per processor patch
    //  procPatchRemoteCC_[patchi][localFaceI] = remote cell centre
    List<pointField> procPatchRemoteCC_;

    //- Whether proc boundary maps have been built
    bool procMapsBuilt_;


    // -- Cached per-cell stencil topology (rebuilt on update()) --

    //- For each seed cell that has been queried:
    //  stencilTopology_[seedCelli] = labelList of size N_,
    //  where entry [a3(i,j,k)] = global cell ID at that stencil position,
    //  or -1 if not found.
    Map<labelList> stencilTopology_;


    // Private Member Functions

    //- Build processor boundary face maps (exchange global IDs and
    //  cell centres across processor patches)
    void buildProcMaps();

    //- Check mesh topology version and rebuild caches if needed.
    //  Called automatically at the start of getZoneField.
    void updateIfNeeded();

    //- Run BFS from seedCelli and return stencil topology:
    //  a labelList of size N_ with global IDs at each stencil position.
    //  Also accumulates remote cell requests into commNeeded.
    labelList buildStencilTopology
    (
        const label seedCelli,
        List<labelHashSet>& commNeeded
    ) const;

    //- Compute relative (i,j,k) offset of targetCC relative to seedCC,
    //  accounting for cyclic wrapping
    Vector<label> relativeIJK
    (
        const point& seedCC,
        const point& targetCC
    ) const;

    //- Check if relative (i,j,k) is within stencil bounds
    bool withinStencil(const Vector<label>& rel) const;

    //- No copy construct
    amrStencil(const amrStencil&) = delete;

    //- No copy assignment
    void operator=(const amrStencil&) = delete;


public:

    // Constructors

    //- Construct from mesh parameters and stencil size
    amrStencil
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
    );

    //- Destructor
    virtual ~amrStencil() = default;


    // Member Functions

    //- Rebuild all cached data (proc maps + stencil topology).
    //  Call after mesh changes (AMR refine/coarsen).
    void update();

    //- True when a stencil topology is cached for this seed cell
    bool hasStencil(const label seedCelli) const
    {
        return stencilTopology_.found(seedCelli);
    }

    //- Populate the stencil for a given seed cell.
    //  Uses cached topology if available, otherwise builds and caches.
    //  phiIJK: pre-collected field map (globalId -> value)
    //  seedCelli: local cell index of the interface cell
    void setStencil
    (
        const Map<scalar>& phiIJK,
        const label seedCelli
    );

    //- Collect field values for all interface cells' stencils.
    //  Uses cached topology to determine which cells are needed,
    //  then exchanges remote values with neighbour processors.
    //  zone: boolean mask for interface cells
    //  phiMap: output map (globalId -> field value)
    //  phi: the volume field to sample
    template<typename Type>
    void getZoneField
    (
        const boolList& zone,
        Map<Type>& phiMap,
        const GeometricField<Type, fvPatchField, volMesh>& phi
    );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "amrStencilTemplates.C"
#endif

#endif

// ************************************************************************* //
