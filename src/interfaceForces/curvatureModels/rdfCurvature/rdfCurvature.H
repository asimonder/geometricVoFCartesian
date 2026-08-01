/*---------------------------------------------------------------------------*\
            Copyright (c) 2026, Asim Onder
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
    Foam::rdfCurvature

Description
    Curvature model based on the Reconstructed Distance Function (RDF).
    K is computed as -div(grad(RDF)/|grad(RDF)|), analogous to divGradAlpha
    but using the smoother signed-distance field from the interface
    reconstruction instead of the raw alpha field. Gives substantially
    more accurate curvature than divGradAlpha while remaining robust on
    non-uniform meshes (unlike the height function).

    Uses the standard OpenFOAM reconstructedDistanceFunction class with
    plicRDF's centre and normal fields.

SourceFiles
    rdfCurvature.C

\*---------------------------------------------------------------------------*/

#ifndef rdfCurvature_H
#define rdfCurvature_H

#include "typeInfo.H"
#include "volFields.H"
#include "dimensionedScalar.H"
#include "autoPtr.H"


#include "curvatureModel.H"

namespace Foam
{

class rdfCurvature
:
    public curvatureModel
{
private:

    //- Disallow copy construct
    rdfCurvature(const rdfCurvature&);

    //- Disallow default bitwise assignment
    void operator=(const rdfCurvature&);

    const dimensionedScalar deltaN_;

    //- RDF field (owned; updated each calculateK)
    //mutable autoPtr<reconstructedDistanceFunction> rdfPtr_;

    //- Band width for RDF construction
    label nBand_;

  const scalar interfaceTol_;
  
protected:

    virtual void calculateK();

    virtual void correctContactAngle
    (
        surfaceVectorField::Boundary& nHat,
        const surfaceVectorField::Boundary& gradAlphaf
    ) const;


public:

    TypeName("rdfCurvature");

    rdfCurvature
    (
        const dictionary& dict,
        const volScalarField& alpha1,
        const volVectorField& U
    );

  // Static utility: compute RDF-based curvature without instantiating
    // a full curvatureModel. Returns a tmp<volScalarField>.
    static tmp<volScalarField> computeKRDF
    (
        const fvMesh& mesh,
        const volScalarField& alpha1
    );

    virtual ~rdfCurvature() = default;
};

  // Static utility: compute RDF-based curvature without instantiating
  // a full curvatureModel. Returns a tmp<volScalarField>.
    static tmp<volScalarField> computeKRDF
    (
        const fvMesh& mesh,
        const volScalarField& alpha1
    );
  
} // End namespace Foam

#endif
