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
    Foam::mlpCurvature

Description


SourceFiles
    mlpCurvature.C

\*---------------------------------------------------------------------------*/

#ifndef mlpCurvature_H
#define mlpCurvature_H

#include "typeInfo.H"
#include "volFields.H"
#include "ijkZone.H"
#include "dimensionedScalar.H"
#include "autoPtr.H"
#include "curvatureModel.H"
#include "globalIndex.H"
#include "reconstructedDistanceFunction.H"
#include "uniformStencil.H"
#include "multilayerPerceptron.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class mlpCurvature Declaration
\*---------------------------------------------------------------------------*/

class mlpCurvature
    : public curvatureModel
{

private:

    //-
  //- Stabilisation for normalisation of the interface normal
  //const dimensionedScalar deltaN_;

  const fvMesh& mesh_;
    // Private Member Functions

  const dimensionedScalar deltaN_;

  multilayerPerceptron mlp_;

  reconstructedDistanceFunction RDF_;
  
  label nMax_;

  boolList boundaryCells_;

  ijkZone ijkMesh_;

  uniformStencil stencil_;
  
  const globalIndex globalNumbering_;

  bool is2D_;

  scalar dl_;

  int iMax_;int jMax_;int kMax_;

  int NInput_;

  //bool bias_;

  bool zonalModel_;

  std::vector<int> indices_;
  
  //- Disallow copy construct
  mlpCurvature(const mlpCurvature&);

  //- Disallow default bitwise assignment
  void operator=(const mlpCurvature&);

  
protected:

  scalar estimateSignK(const List<scalar>& A, const vector& n);
    
  //- Re-calculate the interface curvature
  virtual void calculateK();

  virtual void correctContactAngle
  (
   surfaceVectorField::Boundary& nHat,
   const surfaceVectorField::Boundary& gradAlphaf
   ) const;


public:

    //- Runtime type information
    TypeName("mlpCurvature");


    // Constructors
        //- Construct from components
        mlpCurvature
        (
            const dictionary& dict,
            const volScalarField& alpha1,
            const volVectorField& U
        );

    //- Destructor
    virtual ~mlpCurvature()
    {}

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
