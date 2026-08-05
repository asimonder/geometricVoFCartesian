/*---------------------------------------------------------------------------*\
            Copyright (c) 2026 Asim Onder
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
    Foam::heightFunctionGuarded

Description
    Height-function curvature with guarded columns. Clone of
    heightFunction with four guards on every column of the stencil:

    1. bracket: the column must end pure on both sides
       (alpha < tolPure on the gas end, alpha > 1-tolPure on the
       liquid end), all columns sharing one orientation;
    2. single interface: at most one 0.5-crossing per column;
    3. direction retry: on failure the remaining directions are tried
       in decreasing order of the interface-normal component;
    4. window growth: on failure the window grows symmetrically from
       +-nMax to +-nMaxExtend.

    A cell that passes no guard takes the average curvature of valid
    neighbours in its 3^3 ijk neighbourhood; rdfCurvature remains the
    last resort, as in heightFunction.

    heightFunctionParams: nMax (3), nMaxExtend (8), tolPure (1e-9).

SourceFiles
    heightFunctionGuarded.C

\*---------------------------------------------------------------------------*/

#ifndef heightFunctionGuarded_H
#define heightFunctionGuarded_H

#include "typeInfo.H"
#include "volFields.H"
#include "ijkZone.H"
#include "dimensionedScalar.H"
#include "autoPtr.H"
#include "curvatureModel.H"
#include "globalIndex.H"
#include "zoneDistribute.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class heightFunctionGuarded Declaration
\*---------------------------------------------------------------------------*/

class heightFunctionGuarded
    : public curvatureModel
{

private:

  const fvMesh& mesh_;

  const dimensionedScalar deltaN_;

  ijkZone& ijkMesh_;

  const globalIndex globalNumbering_;

  label nMax_;

  label nMaxExtend_;

  scalar tolPure_;

  bool is2D_;

  scalar dl_;

  scalar interfaceTol_;

  boolList boundaryCells_;

  //- Disallow copy construct
  heightFunctionGuarded(const heightFunctionGuarded&);

  //- Disallow default bitwise assignment
  void operator=(const heightFunctionGuarded&);

  //- Column status of one guarded attempt
  enum tryStatus { OK, IMPURE, MISSING };

  //- Wrapped/reflected lattice index along one axis, replicating the
  //  index rules of heightFunction::calculateHeights
  label wrapX(label b, label o) const;
  label wrapY(label b, label o) const;
  label wrapZ(label b, label o) const;

  //- Guarded column sums over window +-m along dir. On OK fills H
  //  (3 slots in 2D, 9 in 3D) with heights scaled as -sum*dl.
  tryStatus tryHeights
  (
      const Map<scalar>& phiIJK,
      const Vector<label>& ijk,
      const label dir,
      const label m,
      List<scalar>& H
  ) const;

  //- Curvature from the height slots (2D or 3D formula)
  scalar kappaFromHeights(const List<scalar>& H) const;

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
    TypeName("heightFunctionGuarded");


    // Constructors
        //- Construct from components
        heightFunctionGuarded
        (
            const dictionary& dict,
            const volScalarField& alpha1,
            const volVectorField& U
        );

    //- Destructor
    virtual ~heightFunctionGuarded()
    {}

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
