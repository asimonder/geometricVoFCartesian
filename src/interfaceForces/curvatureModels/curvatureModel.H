/*---------------------------------------------------------------------------*\
            Copyright (c) 2021, Asim Onder
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
    Foam::curvatureModel

Description

SourceFiles
    curvatureModel.C
    newCurvatureModel.C

\*---------------------------------------------------------------------------*/

#ifndef curvatureModel_H
#define curvatureModel_H

#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"
#include "dimensionedScalar.H"
#include "autoPtr.H"
#include "Pair.H"

#include "IOdictionary.H"
#include "volFields.H"
#include "surfaceFields.H"



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class curvatureModel Declaration
\*---------------------------------------------------------------------------*/

class curvatureModel
{

protected:

        //- reference to volume fraction field
        const volScalarField& alpha1_;

        //- reference to velcoity
        const volVectorField& U_;

        //- curvature field
        volScalarField K_;
  surfaceScalarField Kf_;


private:

    // Private Member Functions

        //- Disallow copy construct
        curvatureModel(const curvatureModel&);

        //- Disallow default bitwise assignment
        void operator=(const curvatureModel&);



protected:

  //- Re-calculate the interface curvature
  virtual void calculateK(){};

  virtual void correctContactAngle
  (
   surfaceVectorField::Boundary& nHat,
   const surfaceVectorField::Boundary& gradAlphaf
   ) const {};
  
public:

    //- Runtime type information
    TypeName("curvatureModel");

    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            curvatureModel,
            components,
            (
                const dictionary& dict,
                const volScalarField& alpha1,
                const volVectorField& U
            ),
            (dict,alpha1,U)
        );

    // Selectors

        //- Return a reference to the selected phaseChange model
        static autoPtr<curvatureModel> New
        (
            const dictionary& dict,
            const volScalarField& alpha1,
            const volVectorField& U
        );

    // Constructors

        //- Construct from components
        curvatureModel
        (
            const word& type,
            const dictionary& dict,
            const volScalarField& alpha1,
            const volVectorField& U
        );

    //- Destructor
    virtual ~curvatureModel()
    {}

    // Member Functions

  const volScalarField& K() const
  {
    return K_;
  }

  const surfaceScalarField& Kf() const
  {
    return Kf_;
  }

  void correct()
  {
    calculateK();
  }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
