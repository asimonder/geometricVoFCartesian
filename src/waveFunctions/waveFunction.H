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
    Foam::waveFunction

Description

SourceFiles
    waveFunction.C

\*---------------------------------------------------------------------------*/

#ifndef waveFunction_H
#define waveFunction_H

//#include "implicitFunction.H"
#include "autoPtr.H"
#include "dictionary.H"
#include "vector.H"
#include "symmTensor2D.H"
#include "runTimeSelectionTables.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
  //namespace implicitFunctions
  //{

/*---------------------------------------------------------------------------*\
                     Class waveFunction Declaration
\*---------------------------------------------------------------------------*/

class waveFunction
    //    public implicitFunction
{

public:

  //- Runtime type information
  TypeName("waveFunction");

  //- Declare run-time constructor selection table
  declareRunTimeSelectionTable
  (
   autoPtr,
   waveFunction,
   dict,
   (
    const dictionary& dict
    ),
   (dict)
   );

  // Constructors
  waveFunction()=default;

  //- Return a reference to the selected waveFunction
  static autoPtr<waveFunction> New
  (
   const word& waveFunctionType,
   const dictionary& dict
   );

    
  //- Destructor
  virtual ~waveFunction() = default;


  // Member Functions

  virtual scalar value(const vector& p) const
  {
    return GREAT;
  }

  virtual vector grad(const vector& p) const
  {
    return vector::max;
  }

  virtual scalar distanceToSurfaces(const vector& p) const
  {
    return GREAT;
  }

  virtual symmTensor2D hessian(const vector& p) const
  {
    return symmTensor2D::max;
  }
  
  virtual scalar curvature(const vector& p) const
  {
    return GREAT;
  }

  virtual word name() const
  {
    return "NIL";
  }

  virtual scalar getWaterLevel(){return 0.0;};
  virtual void setWaterLevel(scalar H){};
  
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
