/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2025 Asim Onder
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

#include "momentumConsistentFavreMinmod.H"
#include "ITstream.H"
#include "interfaceBand.H"
#include "linear.H"
#include "fvMesh.H"
#include "Pstream.H"
#include "ops.H"

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::surfaceInterpolationScheme<Type>>
Foam::momentumConsistentFavreMinmod<Type>::newBulkScheme
(
    const fvMesh& mesh,
    const surfaceScalarField* faceFluxPtr,
    Istream& is
)
{
    // Copy the remaining tokens into a private stream; the fvSchemes
    // ITstream is cached and rewound between lookups. A leading numeric
    // token (the legacy interfaceThreshold argument) is discarded.
    DynamicList<token> toks;
    while (!is.eof())
    {
        token t(is);
        if (!t.good())
        {
            break;
        }
        if (toks.empty() && t.isNumber())
        {
            continue;   // legacy threshold, ignored
        }
        toks.append(std::move(t));
    }

    if (toks.empty())
    {
        return tmp<surfaceInterpolationScheme<Type>>
        (
            new linear<Type>(mesh)
        );
    }

    ITstream bulkIs(toks);

    if (faceFluxPtr)
    {
        return surfaceInterpolationScheme<Type>::New(mesh, *faceFluxPtr, bulkIs);
    }

    return surfaceInterpolationScheme<Type>::New(mesh, bulkIs);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::surfaceScalarField>
Foam::momentumConsistentFavreMinmod<Type>::weights
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    const fvMesh& mesh = this->mesh();

    // Bulk weights (deep copy; some schemes return non-owning tmps)
    tmp<surfaceScalarField> tWeights
    (
        new surfaceScalarField(tbulkScheme_().weights(vf))
    );
    surfaceScalarField& weights = tWeights.ref();

    const surfaceScalarField& faceMask =
        interfaceBand::New(mesh).faceMask();

    const surfaceScalarField& faceFlux = this->faceFlux_;
    const surfaceScalarField& magSf = mesh.magSf();

    // Statistics counters
    scalar bulkArea = 0.0;
    scalar totalArea = 0.0;

    // Internal faces: the implicit part of the band is pure upwind; the
    // Favre-minmod increment enters explicitly through correction().
    forAll(weights, facei)
    {
        const scalar faceArea = magSf[facei];
        totalArea += faceArea;

        const scalar b = faceMask[facei];
        if (b > SMALL)
        {
            const scalar wUpwind = (faceFlux[facei] > 0) ? 1.0 : 0.0;
            weights[facei] = b*wUpwind + (1.0 - b)*weights[facei];
        }
        bulkArea += (1.0 - b)*faceArea;
    }

    // Boundary faces
    surfaceScalarField::Boundary& bWeights = weights.boundaryFieldRef();
    const surfaceScalarField::Boundary& bFaceMask = faceMask.boundaryField();
    const surfaceScalarField::Boundary& bFaceFlux = faceFlux.boundaryField();
    const surfaceScalarField::Boundary& bMagSf = magSf.boundaryField();

    forAll(bWeights, patchi)
    {
        scalarField& pWeights = bWeights[patchi];
        const scalarField& pFaceMask = bFaceMask[patchi];
        const scalarField& pFaceFlux = bFaceFlux[patchi];
        const scalarField& pMagSf = bMagSf[patchi];

        forAll(pWeights, facei)
        {
            const scalar faceArea = pMagSf[facei];
            totalArea += faceArea;

            const scalar b = pFaceMask[facei];
            if (b > SMALL)
            {
                const scalar wUpwind = (pFaceFlux[facei] > 0) ? 1.0 : 0.0;
                pWeights[facei] = b*wUpwind + (1.0 - b)*pWeights[facei];
            }
            bulkArea += (1.0 - b)*faceArea;
        }
    }

    // Parallel reduction of the statistics
    reduce(bulkArea, sumOp<scalar>());
    reduce(totalArea, sumOp<scalar>());

    if (Pstream::master())
    {
        scalar bulkFraction = 0.0;
        if (totalArea > VSMALL)
        {
            bulkFraction = bulkArea / totalArea;
        }

        Info<< "MomentumFavreMinmod: Bulk Face Area Fraction = "
            << bulkFraction * 100.0 << " %" << endl;
    }

    return tWeights;
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh>>
Foam::momentumConsistentFavreMinmod<Type>::correction
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    const fvMesh& mesh = this->mesh();

    const surfaceScalarField& faceMask =
        interfaceBand::New(mesh).faceMask();
    const surfaceScalarField& faceFlux = this->faceFlux_;

    if (!mesh.template foundObject<volScalarField>("rho"))
    {
        FatalErrorInFunction
            << "momentumConsistentFavreMinmod requires the registered"
            << " density field rho" << nl
            << exit(FatalError);
    }
    const volScalarField& rho =
        mesh.template lookupObject<volScalarField>("rho");

    // Favre-averaged minmod face value, Kuhn et al. (2023) Eq. 35 with
    // one SHARED limiter computed on rho: interp_w(rho*vf)/interp_w(rho)
    // with the scalar-minmod weights w of rho applied to both fields.
    // Per-field limiters are not usable here: the stock vector limiters
    // evaluate their ratio on magSqr of the field, so the rho*U limiter
    // would see rho^2 while the rho limiter sees rho, and the ratio then
    // fails to cancel for uniform flow (measured Galilean violation).
    // With one shared weight the ratio is exact for uniform flow, the
    // face value stays a convex combination of the two cell values, and
    // the density jump is carried by the mass weighting (the continuous-
    // velocity argument of Arrufat et al. 2021, Eq. 48).
    DynamicList<token> tokB;
    tokB.append(token(word("Minmod")));
    ITstream isB(tokB);
    tmp<surfaceInterpolationScheme<scalar>> tMinmodS
    (
        surfaceInterpolationScheme<scalar>::New(mesh, faceFlux, isB)
    );
    tmp<surfaceScalarField> tw(tMinmodS().weights(rho));
    const surfaceScalarField& w = tw();

    tmp<GeometricField<Type, fvPatchField, volMesh>> trhoVf(rho*vf);
    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> tNum
    (
        surfaceInterpolationScheme<Type>::interpolate(trhoVf(), w)
    );
    trhoVf.clear();

    tmp<surfaceScalarField> tDen
    (
        surfaceInterpolationScheme<scalar>::interpolate(rho, w)
    );

    // Pure upwind face value of vf (the implicit part in the band)
    upwind<Type> upwindScheme(mesh, faceFlux);
    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> tUpwindVf
    (
        upwindScheme.interpolate(vf)
    );

    // Band correction towards the Favre-minmod value
    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> tcorr
    (
        faceMask*(tNum()/tDen() - tUpwindVf())
    );

    // Bulk correction outside the band
    if (tbulkScheme_().corrected())
    {
        tcorr.ref() +=
            (1.0 - faceMask)*tbulkScheme_().correction(vf);
    }

    return tcorr;
}


namespace Foam
{
    makeSurfaceInterpolationScheme(momentumConsistentFavreMinmod)
}

// ************************************************************************* //
