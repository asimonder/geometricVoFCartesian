/*---------------------------------------------------------------------------*\
            Copyright (c) 2021 Asim Onder
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
    Foam::heightFunction

Description


SourceFiles
    heightFunction.C

\*---------------------------------------------------------------------------*/

#ifndef heightFunction_H
#define heightFunction_H

#include "typeInfo.H"
#include "volFields.H"
#include "ijkZone.H"
#include "dimensionedScalar.H"
#include "autoPtr.H"
#include "curvatureModel.H"
#include "globalIndex.H"
//#include "reconstructedDistanceFunction.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class heightFunction Declaration
\*---------------------------------------------------------------------------*/

class heightFunction
    : public curvatureModel
{

private:

    //-
  //- Stabilisation for normalisation of the interface normal
  //const dimensionedScalar deltaN_;

  const fvMesh& mesh_;
    // Private Member Functions

  const dimensionedScalar deltaN_;

  ijkZone& ijkMesh_;

  const globalIndex globalNumbering_;

  label nMax_;

  bool is2D_;

  scalar dl_;

  scalar interfaceTol_;

  boolList boundaryCells_;

  
  //- Disallow copy construct
  heightFunction(const heightFunction&);

  //- Disallow default bitwise assignment
  void operator=(const heightFunction&);

  List<scalar> calculateHeights(const Map<scalar>& phiIJK,const Vector<label>& ijk,label dir);
  
protected:

  //- Re-calculate the interface curvature
  virtual void calculateK();

  virtual void correctContactAngle
  (
   surfaceVectorField::Boundary& nHat,
   const surfaceVectorField::Boundary& gradAlphaf
   ) const;


public:

    //- Runtime type information
    TypeName("heightFunction");


    // Constructors
        //- Construct from components
        heightFunction
        (
            const dictionary& dict,
            const volScalarField& alpha1,
            const volVectorField& U
        );

    //- Destructor
    virtual ~heightFunction()
    {}

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
