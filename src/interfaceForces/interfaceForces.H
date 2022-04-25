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
    Foam::interfaceForces

Description

SourceFiles
    interfaceForces.C

\*---------------------------------------------------------------------------*/

#ifndef interfaceForces_H
#define interfaceForces_H

#include "IOdictionary.H"
#include "surfaceTensionModel.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "curvatureModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class interfaceForces Declaration
\*---------------------------------------------------------------------------*/

class interfaceForces
{
  // Private data

  //- Keep a reference to the transportProperties dictionary
  const dictionary& transportPropertiesDict_;

  //- Surface tension
  autoPtr<surfaceTensionModel> sigmaPtr_;


  autoPtr<curvatureModel> curvature_;

  //- Stabilisation for normalisation of the interface normal
  //const dimensionedScalar deltaN_;
  
  const volScalarField& alpha1_;
  const volVectorField& U_;
  //surfaceScalarField nHatf_;
  //volScalarField K_;
  

  // Private Member Functions

  //- No copy construct
  interfaceForces(const interfaceForces&) = delete;
  
  //- No copy assignment
  void operator=(const interfaceForces&) = delete;

  //- Correction for the boundary condition on the unit normal nHat on
  //  walls to produce the correct contact dynamic angle
  //  calculated from the component of U parallel to the wall
  /*void correctContactAngle
        (
	 surfaceVectorField::Boundary& nHat,
	 const surfaceVectorField::Boundary& gradAlphaf
	 ) const;

  
  //protected:
  //- Re-calculate the interface curvature
  void calculateK();*/


public:

    // Constructors

        //- Construct from volume fraction field gamma and IOdictionary
  interfaceForces
  (
   const volScalarField& alpha1,
   const volVectorField& U,
   const IOdictionary&
   );
  
  
  // Member Functions
  
  /*scalar cAlpha() const
  {
    return cAlpha_;
  }

  const dimensionedScalar& deltaN() const
  {
    return deltaN_;
  }

  const surfaceScalarField& nHatf() const
  {
    return curvature_->nHatf_;
    }*/

  //tmp<volScalarField> sigmaK() const;
  tmp<volScalarField> sigma() const;

  tmp<volScalarField> K() const;

  tmp<surfaceScalarField> Kf() const;
  
  tmp<surfaceScalarField> surfaceTensionForce() const;
  //tmp<volVectorField> surfaceTensionForce() const;

  //- Indicator of the proximity of the interface
  //  Field values are 1 near and 0 away for the interface.
  tmp<volScalarField> nearInterface() const;

  //- Read transportProperties dictionary
  bool read();

  void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
