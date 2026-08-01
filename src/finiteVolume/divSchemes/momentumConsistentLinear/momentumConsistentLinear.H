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

Class
    Foam::momentumConsistentLinear

Description
    Interface-band blended convection interpolation for mass-momentum
    consistency with geometric VOF (isoAdvector): forces upwind weights on
    faces inside the interface band (provided by the interfaceBand mesh
    object, controlled by interfaceThreshold/nBufferLayers in fvSolution),
    and uses an arbitrary bulk scheme read from the fvSchemes line
    elsewhere.

    Usage:
    \verbatim
    divSchemes
    {
        div(rhoPhi,U)  Gauss momentumConsistentLinear limitedLinearV 1;
    }
    \endverbatim

    With no bulk scheme given, linear is used (the original behaviour).
    A leading numeric token (the legacy interfaceThreshold argument) is
    accepted and ignored; the threshold now lives in fvSolution.

\*---------------------------------------------------------------------------*/

#include "upwind.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifndef momentumConsistentLinear_H
#define momentumConsistentLinear_H

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class momentumConsistentLinear Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class momentumConsistentLinear
:
    public upwind<Type>
{
    // Private Data

        //- Bulk interpolation scheme (used outside the interface band)
        tmp<surfaceInterpolationScheme<Type>> tbulkScheme_;


    // Private Member Functions

        //- Construct the bulk scheme from the remaining stream tokens;
        //  swallows a legacy numeric threshold argument; defaults to linear
        static tmp<surfaceInterpolationScheme<Type>> newBulkScheme
        (
            const fvMesh& mesh,
            const surfaceScalarField* faceFluxPtr,
            Istream& is
        );

        //- No copy construct
        momentumConsistentLinear(const momentumConsistentLinear&) = delete;

        //- No copy assignment
        void operator=(const momentumConsistentLinear&) = delete;


public:

    //- Runtime type information
    TypeName("momentumConsistentLinear");


    // Constructors

        //- Construct from mesh, faceFlux and Istream
        momentumConsistentLinear
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            Istream& schemeData
        )
        :
            upwind<Type>(mesh, faceFlux, schemeData),
            tbulkScheme_(newBulkScheme(mesh, &faceFlux, schemeData))
        {}

        //- Construct from mesh and Istream
        momentumConsistentLinear
        (
            const fvMesh& mesh,
            Istream& schemeData
        )
        :
            upwind<Type>(mesh, schemeData),
            tbulkScheme_(newBulkScheme(mesh, nullptr, schemeData))
        {}


    // Member Functions

        //- Blended interpolation weights
        virtual tmp<surfaceScalarField> weights
        (
            const GeometricField<Type, fvPatchField, volMesh>& vf
        ) const;

        //- True if the bulk scheme carries an explicit correction
        //  (e.g. linearUpwind); the correction is masked in the band
        virtual bool corrected() const
        {
            return tbulkScheme_().corrected();
        }

        //- Bulk-scheme correction, zeroed on interface-band faces so the
        //  band stays pure upwind while the bulk keeps its full order
        virtual tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
        correction
        (
            const GeometricField<Type, fvPatchField, volMesh>& vf
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
