/*---------------------------------------------------------------------------*\
            Copyright (c) 2021, Asim Onder
-------------------------------------------------------------------------------
License
    This file is part of the VoFLibrary source code library, which is an
	unofficial extension to OpenFOAM.
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
    Foam::divGradAlpha

Description


SourceFiles
    divGradAlpha.C

\*---------------------------------------------------------------------------*/

#ifndef divGradAlpha_H
#define divGradAlpha_H

#include "typeInfo.H"
#include "volFields.H"
#include "dimensionedScalar.H"
#include "autoPtr.H"




#include "curvatureModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class divGradAlpha Declaration
\*---------------------------------------------------------------------------*/

class divGradAlpha
    : public curvatureModel
{

private:

    // Private Member Functions

  //- Disallow copy construct
  divGradAlpha(const divGradAlpha&);
  
  //- Disallow default bitwise assignment
  void operator=(const divGradAlpha&);

  const dimensionedScalar deltaN_;

protected:

        //- Calculate the interface curvature
  virtual void calculateK();
  
  virtual void correctContactAngle
  (
   surfaceVectorField::Boundary& nHat,
   const surfaceVectorField::Boundary& gradAlphaf
   ) const;



public:

    //- Runtime type information
    TypeName("divGradAlpha");

    // Constructors

  //- Construct from components
  divGradAlpha
  (
   const dictionary& dict,
   const volScalarField& alpha1,
   const volVectorField& U
   );


  //- Destructor
  virtual ~divGradAlpha()
  {}


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
