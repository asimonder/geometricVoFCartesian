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
#include "fvMesh.H"
#include "fvcSnGrad.H"
// for parallel reduction
#include "Pstream.H"
#include "ops.H" 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::surfaceScalarField>
Foam::momentumConsistentLinear<Type>::weights
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    const fvMesh& mesh = this->mesh();

    tmp<surfaceScalarField> tWeights
    (
        new surfaceScalarField
        (
            IOobject
            (
                "momentumConsistentLinear::weights",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            mesh,
            dimless
        )
    );

    surfaceScalarField& weights = tWeights.ref();

    tmp<surfaceScalarField> tsnGradAlpha = fvc::snGrad(alpha_);
    const surfaceScalarField& snGradAlpha = tsnGradAlpha();

    const surfaceScalarField& faceFlux = this->faceFlux_;
    const surfaceScalarField& meshWeights = mesh.weights();
    
    // Get face areas for accurate statistics
    const surfaceScalarField& magSf = mesh.magSf();

    // Statistics counters
    scalar linearArea = 0.0;
    scalar totalArea = 0.0;

    // Internal faces
    forAll(weights, facei)
    {
        scalar faceArea = magSf[facei];
        totalArea += faceArea;

        if (mag(snGradAlpha[facei]) > interfaceThreshold_)
        {
            // Interface: use upwind
            weights[facei] = (faceFlux[facei] > 0) ? 1.0 : 0.0;
        }
        else
        {
            // Bulk: use linear
            weights[facei] = meshWeights[facei];
            linearArea += faceArea; // Count this area as Linear/Bulk
        }
    }

    // Boundary faces
    surfaceScalarField::Boundary& bWeights = weights.boundaryFieldRef();
    const surfaceScalarField::Boundary& bSnGradAlpha = snGradAlpha.boundaryField();
    const surfaceScalarField::Boundary& bFaceFlux = faceFlux.boundaryField();
    const surfaceScalarField::Boundary& bMeshWeights = meshWeights.boundaryField();
    const surfaceScalarField::Boundary& bMagSf = magSf.boundaryField();

    forAll(bWeights, patchi)
    {
        scalarField& pWeights = bWeights[patchi];
        const scalarField& pSnGradAlpha = bSnGradAlpha[patchi];
        const scalarField& pFaceFlux = bFaceFlux[patchi];
        const scalarField& pMeshWeights = bMeshWeights[patchi];
        const scalarField& pMagSf = bMagSf[patchi];

        forAll(pWeights, facei)
        {
            scalar faceArea = pMagSf[facei];
            totalArea += faceArea;

            if (mag(pSnGradAlpha[facei]) > interfaceThreshold_)
            {
                pWeights[facei] = (pFaceFlux[facei] > 0) ? 1.0 : 0.0;
            }
            else
            {
                pWeights[facei] = pMeshWeights[facei];
                linearArea += faceArea; // Count as Linear
            }
        }
    }

    // Parallel Reduction
    // This ensures the stats are summed across all processors
    reduce(linearArea, sumOp<scalar>());
    reduce(totalArea, sumOp<scalar>());

    // Print Statistics (Only on Master to avoid log spam)
    if (Pstream::master())
    {
        scalar bulkFraction = 0.0;
        if (totalArea > VSMALL)
        {
            bulkFraction = linearArea / totalArea;
        }
        
        // Print the fraction of the Face Area that is using Linear interpolation
        Info<< "MomentumLinear: Bulk (Linear) Face Area Fraction = " 
            << bulkFraction * 100.0 << " %" << endl;
    }

    return tWeights;
}

namespace Foam
{
    makeSurfaceInterpolationScheme(momentumConsistentLinear)
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //







