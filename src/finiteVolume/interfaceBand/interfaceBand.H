/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2026 Asim Onder
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
    Foam::interfaceBand

Description
    Mesh-registered singleton (wallDist/zoneDistribute pattern) providing
    the interface band for mass-momentum consistent blended schemes:
    band = mixed cells (threshold < alpha < 1-threshold) plus nBufferLayers
    face-connected rings around them.

    The rings are propagated with FaceCellWave<topoDistanceData<label>>,
    which is correct across processor and cyclic boundaries; the wave cost
    scales with the band, not the mesh. Rebuilt lazily when the time index
    changes, so the first user in each time step (momentum assembly, which
    runs after the alpha advection) sees the post-advection band. Deleted
    and rebuilt automatically on mesh topology change (AMR).

    Controls, read from the interfaceBand sub-dictionary of fvSolution and
    re-read on every rebuild (so runTimeModifiable picks up changes):
    \verbatim
    interfaceBand
    {
        nLayers     3;      // hard rings beyond mixed cells (band = 1)
        nRampLayers 0;      // graded rings beyond the core (cosine decay)
        threshold   1e-6;   // mixed-cell detection on alpha
        alpha       alpha.water;  // optional; auto-detected if absent
    }
    \endverbatim
    Legacy flat fvSolution keys (interfaceThreshold, nBufferLayers,
    interfaceBandAlpha) are accepted as fallback; the sub-dictionary wins.

    The cell band is registered (and written) as the volScalarField
    "interfaceBandCells" for inspection.

SourceFiles
    interfaceBand.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_interfaceBand_H
#define Foam_interfaceBand_H

#include "MeshObject.H"
#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class interfaceBand Declaration
\*---------------------------------------------------------------------------*/

class interfaceBand
:
    public MeshObject<fvMesh, TopologicalMeshObject, interfaceBand>
{
    // Private Data

        //- Time index of the last rebuild (-1: never built)
        label timeIndex_;

        //- Number of face-connected buffer rings beyond mixed cells
        label nLayers_;

        //- Number of graded transition rings beyond the hard core
        label nRampLayers_;

        //- Mixed-cell detection threshold on alpha
        scalar threshold_;

        //- Name of the volume-fraction field
        word alphaName_;

        //- Cell marker: 1 in mixed cells + nLayers_ rings, 0 in bulk
        volScalarField band_;

        //- Face mask: 1 where either adjacent cell is in the band
        surfaceScalarField faceMask_;


    // Private Member Functions

        //- Auto-detect the registered volume-fraction field name
        static word findAlphaName(const fvMesh& mesh);

        //- (Re-)read controls from fvSolution; logs when values change
        void readControls();

        //- Rebuild band and face mask from the current alpha
        void rebuild();

        //- Rebuild if the time index changed since the last build
        void updateIfNeeded();

        //- No copy construct
        interfaceBand(const interfaceBand&) = delete;

        //- No copy assignment
        void operator=(const interfaceBand&) = delete;


public:

    //- Runtime type information
    TypeName("interfaceBand");


    // Constructors

        //- Construct from mesh; reads controls from fvSolution
        explicit interfaceBand(const fvMesh&);


    //- Selector
    static interfaceBand& New(const fvMesh&);


    //- Destructor
    virtual ~interfaceBand() = default;


    // Member Functions

        //- Cell band marker, updated to the current time step
        const volScalarField& band();

        //- Face mask, updated to the current time step
        const surfaceScalarField& faceMask();

        //- Number of hard buffer rings
        label nLayers() const noexcept
        {
            return nLayers_;
        }

        //- Number of graded transition rings
        label nRampLayers() const noexcept
        {
            return nRampLayers_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
