/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2026 Asim Onder
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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

\*---------------------------------------------------------------------------*/

#include "midpointGaussConvectionScheme.H"
#include "fvcSurfaceIntegrate.H"
#include "fvMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fv
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
const surfaceInterpolationScheme<Type>&
midpointGaussConvectionScheme<Type>::interpScheme() const
{
    return tinterpScheme_();
}


template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
midpointGaussConvectionScheme<Type>::interpolate
(
    const surfaceScalarField&,
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    // Trapezoidal average of the transported field.
    // Each level uses the scheme's own weights; for weight-only schemes
    // whose weights depend on flux sign and band membership (upwind,
    // momentumConsistentLinear) the two weight fields are identical.
    return
        0.5*
        (
            tinterpScheme_().interpolate(vf)
          + tinterpScheme_().interpolate(vf.oldTime())
        );
}


template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
midpointGaussConvectionScheme<Type>::flux
(
    const surfaceScalarField& faceFlux,
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    return faceFlux*interpolate(faceFlux, vf);
}


template<class Type>
tmp<fvMatrix<Type>>
midpointGaussConvectionScheme<Type>::fvmDiv
(
    const surfaceScalarField& faceFlux,
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    // Single weight field at the current level, shared by both halves.
    tmp<surfaceScalarField> tweights = tinterpScheme_().weights(vf);
    const surfaceScalarField& weights = tweights();

    tmp<fvMatrix<Type>> tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            faceFlux.dimensions()*vf.dimensions()
        )
    );
    fvMatrix<Type>& fvm = tfvm.ref();

    // Implicit half at n+1: stock Gauss coefficients scaled by 0.5
    fvm.lower() = -0.5*weights.primitiveField()*faceFlux.primitiveField();
    fvm.upper() = fvm.lower() + 0.5*faceFlux.primitiveField();
    fvm.negSumDiag();

    forAll(vf.boundaryField(), patchi)
    {
        const fvPatchField<Type>& psf = vf.boundaryField()[patchi];
        const fvsPatchScalarField& patchFlux = faceFlux.boundaryField()[patchi];
        const fvsPatchScalarField& pw = weights.boundaryField()[patchi];

        fvm.internalCoeffs()[patchi] =
            0.5*patchFlux*psf.valueInternalCoeffs(pw);
        fvm.boundaryCoeffs()[patchi] =
           -0.5*patchFlux*psf.valueBoundaryCoeffs(pw);
    }

    if (tinterpScheme_().corrected())
    {
        // Half correction at each level. Dormant for weight-only schemes
        // (upwind-derived blends report corrected() == false).
        fvm +=
            0.5*fvc::surfaceIntegrate
            (
                faceFlux*tinterpScheme_().correction(vf)
            );
        fvm +=
            0.5*fvc::surfaceIntegrate
            (
                faceFlux*tinterpScheme_().correction(vf.oldTime())
            );
    }

    // Explicit half from oldTime() with the SAME weights. Same source
    // sign convention as the stock corrected() branch. This call clears
    // tweights, so it must remain the last use of the weight field.
    fvm +=
        0.5*fvc::surfaceIntegrate
        (
            faceFlux
           *surfaceInterpolationScheme<Type>::interpolate
            (
                vf.oldTime(),
                tweights
            )
        );

    return tfvm;
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
midpointGaussConvectionScheme<Type>::fvcDiv
(
    const surfaceScalarField& faceFlux,
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    tmp<GeometricField<Type, fvPatchField, volMesh>> tConvection
    (
        fvc::surfaceIntegrate(flux(faceFlux, vf))
    );

    tConvection.ref().rename
    (
        "convection(" + faceFlux.name() + ',' + vf.name() + ')'
    );

    return tConvection;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

} // End namespace Foam

// ************************************************************************* //
