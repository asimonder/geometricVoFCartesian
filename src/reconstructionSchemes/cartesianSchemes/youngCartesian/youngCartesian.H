/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 Asim Onder
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
    Foam::reconstruction::youngCartesian

Description


SourceFiles
    youngCartesian.C

\*---------------------------------------------------------------------------*/

#ifndef youngCartesian_H
#define youngCartesian_H

#include "typeInfo.H"
#include "reconstructionSchemes.H"
#include "volFields.H"
#include "dimensionedScalar.H"
#include "autoPtr.H"
#include "surfaceIteratorPLIC.H"
#include "ijkZone.H"
#include "uniformStencil.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace reconstruction
{

/*---------------------------------------------------------------------------*\
                          Class youngCartesian Declaration
\*---------------------------------------------------------------------------*/

class youngCartesian
:
    public reconstructionSchemes
{
    // Private Data

  //- Reference to mesh
  const fvMesh& mesh_;

  //- Interpolation object from cell centres to points
  DynamicField<vector> interfaceNormal_;
  
  ijkZone ijkMesh_;

  uniformStencil stencil_;

  boolList boundaryCells_;

  const globalIndex globalNumbering_;
  
  // Switches and tolerances. Tolerances need to go into toleranceSwitches

  //- Tolerance for search of isoFace giving specified VOF value
  scalar isoFaceTol_;

  //- Tolerance for marking of surface cells:
  //  Those with surfCellTol_ < alpha1 < 1 - surfCellTol_
  scalar surfCellTol_;
  
  //- SurfaceIterator finds the plane centre for specified VOF value
  surfaceIteratorPLIC sIterPLIC_;

  //- Compute gradient at the surfaces
  void gradSurf(const volScalarField& phi);

  // Private Member Functions

  //- No copy construct
  youngCartesian(const youngCartesian&) = delete;
  
  //- No bitwise assignment
  void operator=(const youngCartesian&) = delete;

  vector calcYoungNormal(const scalarList& A);
  //vector calcCCDNormal(const scalarList& A);


public:

    //- Runtime type information
    TypeName("youngCartesian");

    //- Construct from components
    youngCartesian
    (
        volScalarField& alpha1,
        const surfaceScalarField& phi,
        const volVectorField& U,
        const dictionary& dict
    );

    //- Destructor
    virtual ~youngCartesian() = default;


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
