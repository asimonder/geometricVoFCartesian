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

#include "momentumConsistentLinear.H"
#include "ITstream.H"
#include "interfaceBand.H"
#include "linear.H"
#include "fvMesh.H"
// for parallel reduction
#include "Pstream.H"
#include "ops.H"

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::surfaceInterpolationScheme<Type>>
Foam::momentumConsistentLinear<Type>::newBulkScheme
(
    const fvMesh& mesh,
    const surfaceScalarField* faceFluxPtr,
    Istream& is
)
{
    // Copy the remaining tokens into a private stream. The fvSchemes
    // ITstream is cached and rewound between lookups, so it must not be
    // left holding putBack state; and its eof flag hides a putBack token
    // from the scheme selector. A leading numeric token (the legacy
    // interfaceThreshold argument) is discarded.
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
        // No bulk scheme given: default to linear
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
Foam::momentumConsistentLinear<Type>::weights
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

    // Get face areas for accurate statistics
    const surfaceScalarField& magSf = mesh.magSf();

    // Statistics counters
    scalar bulkArea = 0.0;
    scalar totalArea = 0.0;

    // Internal faces. The face factor b is 1 in the core (pure upwind),
    // fractional in the ramp rings, 0 in the bulk.
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

    // Parallel Reduction
    // This ensures the stats are summed across all processors
    reduce(bulkArea, sumOp<scalar>());
    reduce(totalArea, sumOp<scalar>());

    // Print Statistics (Only on Master to avoid log spam)
    if (Pstream::master())
    {
        scalar bulkFraction = 0.0;
        if (totalArea > VSMALL)
        {
            bulkFraction = bulkArea / totalArea;
        }

        // Fraction of the face area using the bulk scheme
        Info<< "MomentumLinear: Bulk Face Area Fraction = "
            << bulkFraction * 100.0 << " %" << endl;
    }

    return tWeights;
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh>>
Foam::momentumConsistentLinear<Type>::correction
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    const fvMesh& mesh = this->mesh();

    // Deep copy of the bulk correction
    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> tcorr
    (
        new GeometricField<Type, fvsPatchField, surfaceMesh>
        (
            tbulkScheme_().correction(vf)
        )
    );

    // Zero the correction on band faces: band faces are pure upwind
    tcorr.ref() *= (scalar(1) - interfaceBand::New(mesh).faceMask());

    return tcorr;
}


namespace Foam
{
    makeSurfaceInterpolationScheme(momentumConsistentLinear)
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
