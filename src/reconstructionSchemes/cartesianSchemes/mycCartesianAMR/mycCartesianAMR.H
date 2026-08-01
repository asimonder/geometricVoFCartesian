/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2026 Asim Onder
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of the unofficial extension library geometricVoFCartesian.

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
    Foam::reconstruction::mycCartesianAMR

Description
    MYC (mixed Youngs-centered) interface reconstruction on adaptively
    refined Cartesian meshes. Identical normal selection to mycCartesian,
    with the stencil collected through amrStencil (BFS neighbour
    discovery with geometric ijk offsets) instead of the global ijkZone
    lookup, following the heightFunctionAMR pattern. Valid because the
    refinement band keeps every interface cell and its 27-cell cube at
    the finest uniform level.

    Cells whose cube is incomplete fall back to the Youngs normal; cells
    at non-cyclic physical boundaries are snapped as in mycCartesian.

SourceFiles
    mycCartesianAMR.C

\*---------------------------------------------------------------------------*/

#ifndef mycCartesianAMR_H
#define mycCartesianAMR_H

#include "typeInfo.H"
#include "reconstructionSchemes.H"
#include "volFields.H"
#include "dimensionedScalar.H"
#include "autoPtr.H"
#include "surfaceIteratorPLIC.H"
#include "amrStencil.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace reconstruction
{

/*---------------------------------------------------------------------------*\
                      Class mycCartesianAMR Declaration
\*---------------------------------------------------------------------------*/

class mycCartesianAMR
:
    public reconstructionSchemes
{
    // Private Data

        //- Reference to mesh
        const fvMesh& mesh_;

        //- Interface normal per interface cell
        DynamicField<vector> interfaceNormal_;

        //- Tolerance for search of isoFace giving specified VOF value
        scalar isoFaceTol_;

        //- Tolerance for marking of surface cells
        scalar surfCellTol_;

        //- SurfaceIterator finds the plane centre for specified VOF value
        surfaceIteratorPLIC sIterPLIC_;

        //- Empty (2D) direction flags
        Vector<bool> isEmpty_;

        //- Domain thickness in the empty direction
        scalar emptyDirThickness_;

        //- AMR-capable stencil (BFS discovery, cached topology)
        autoPtr<amrStencil> amrStencilPtr_;

        //- Cells within one layer of a non-cyclic physical boundary
        boolList boundaryCells_;

        //- Cross-check hash slots against findCell (diagnostic)
        bool debugCheck_;


    // Private Member Functions

        //- Compute interface normals from the 27-cell stencils
        void gradSurf(const volScalarField& phi);

        //- Mark cells within nLayers of non-cyclic physical boundaries
        void markBoundaryCells(label nLayers);

        //- No copy construct
        mycCartesianAMR(const mycCartesianAMR&) = delete;

        //- No bitwise assignment
        void operator=(const mycCartesianAMR&) = delete;


public:

    //- Runtime type information
    TypeName("mycCartesianAMR");

    //- Construct from components
    mycCartesianAMR
    (
        volScalarField& alpha1,
        const surfaceScalarField& phi,
        const volVectorField& U,
        const dictionary& dict
    );

    //- Destructor
    virtual ~mycCartesianAMR() = default;


    // Member Functions

        //- Reconstruct interface
        virtual void reconstruct(bool forceUpdate = true);

        //- map VoF Field in case of refinement
        virtual void mapAlphaField() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace reconstruction
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
