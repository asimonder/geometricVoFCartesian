/*---------------------------------------------------------------------------*\
            Copyright (c) 2021-2026 Asim Onder
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
    Foam::heightFunctionAMR

Description
    Height function curvature using amrStencil for field collection.
    Uses direct position-to-alpha hash for O(1) stencil lookups.
    Direction-aware: 63 lookups per cell (3x3x7), not 343 (7x7x7).
    Targeted face interpolation: O(nInterface) not O(nAllFaces).

    Fallback methods for incomplete stencils (dict selectable):
      none         - K=0 at fallback cells
      rdf          - full computeKRDF (accurate, slow at scale)
      divGradAlpha - div(grad(alpha)) (fast, low accuracy)
      interpolate  - average from face-neighbor valid K cells

SourceFiles
    heightFunctionAMR.C

\*---------------------------------------------------------------------------*/

#ifndef heightFunctionAMR_H
#define heightFunctionAMR_H

#include "typeInfo.H"
#include "volFields.H"
#include "dimensionedScalar.H"
#include "autoPtr.H"
#include "curvatureModel.H"
#include "amrStencil.H"
#include "zoneDistribute.H"

namespace Foam
{

class heightFunctionAMR
    : public curvatureModel
{

private:

  const fvMesh& mesh_;

  const dimensionedScalar deltaN_;

  Vector<bool> isEmpty_;

  scalar emptyDirThickness_;

  label nMax_;

  bool is2D_;

  scalar interfaceTol_;

  boolList boundaryCells_;

  autoPtr<amrStencil> amrStencilPtr_;

  Vector<bool> isCyclic_;

  scalar Lx_;
  scalar Ly_;
  scalar Lz_;

  //- Fallback method: "none", "rdf", "divGradAlpha", "interpolate"
  word fallbackMethod_;

  heightFunctionAMR(const heightFunctionAMR&);

  void operator=(const heightFunctionAMR&);

  void markBoundaryCells(label nLayers);

protected:

  virtual void calculateK();

  virtual void correctContactAngle
  (
   surfaceVectorField::Boundary& nHat,
   const surfaceVectorField::Boundary& gradAlphaf
  ) const;

public:

    TypeName("heightFunctionAMR");

    heightFunctionAMR
    (
        const dictionary& dict,
        const volScalarField& alpha1,
        const volVectorField& U
    );

    virtual ~heightFunctionAMR()
    {}
};

} // End namespace Foam

#endif
