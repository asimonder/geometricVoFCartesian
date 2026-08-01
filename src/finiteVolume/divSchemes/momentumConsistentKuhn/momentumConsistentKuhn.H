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
    Foam::momentumConsistentKuhn

Description
    Interface-band convection interpolation after Kuhn, Deskos, Sprague
    (2023), Computers and Fluids 252:105770, Eqs. 31, 32, 34, 35.

    Since 2026-07-12 the Delta-t tracing terms are INCLUDED (his
    instruction, reversing the 2026-07 spatial-only spec after the
    dt*t secular crest-growth diagnosis): the upstream offset of
    Eq. 31 is shortened by u_adv*dt/2 using the registered volumetric
    flux phi (the alpha-advecting flux, so face values are centred on
    the swept space-time volumes), and the Eq. 34 transverse half-step
    terms are applied to numerator and denominator alike. Deviation
    from the paper: transverse gradients use the same minmod slopes as
    Eq. 32 and a cell-centred fvc::reconstruct(phi) advecting velocity
    instead of their PPM/Almgren machinery (more dissipative, stable).
    A non-positive traced denominator falls back to pure upwind.

    Band faces (interfaceBand mesh object; Kuhn's flag equals
    nLayers 1) are implicit pure upwind plus an explicit deferred
    correction towards the Favre value

        u_f = rec(rho U)_f / rec(rho)_f,

    where rec() is the UPSTREAM-cell reconstruction: upstream value
    plus the minmod of the two one-sided differences of the upstream
    cell along the face direction, times the face offset (Kuhn
    Eqs. 31-32, per field, per component). Both numerator and
    denominator are anchored to the upstream cell, so the face never
    prices the flux with the downwind velocity and the flux cannot
    remove more momentum than the upstream cell holds (Kuhn Fig. 3).
    Minmod is positively homogeneous, so uniform flow cancels the
    correction exactly at any density ratio.

    Slopes need the second-upstream neighbour: computed on axis-aligned
    faces with single plus/minus neighbours per direction; anywhere
    that structure is absent (non-Cartesian faces, refinement edges,
    walls) the slope falls back to zero, which is pure upwind.

    Usage:
    \verbatim
    divSchemes
    {
        div(rhoPhi,U)  Gauss momentumConsistentKuhn linear;
    }
    \endverbatim

    With no bulk scheme given, linear is used. A leading numeric token
    (legacy interfaceThreshold argument) is accepted and ignored.

\*---------------------------------------------------------------------------*/

#include "upwind.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifndef momentumConsistentKuhn_H
#define momentumConsistentKuhn_H

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class momentumConsistentKuhn Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class momentumConsistentKuhn
:
    public upwind<Type>
{
    // Private Data

        //- Bulk interpolation scheme (used outside the interface band)
        tmp<surfaceInterpolationScheme<Type>> tbulkScheme_;


    // Private Member Functions

        //- Construct the bulk scheme from the remaining stream tokens
        static tmp<surfaceInterpolationScheme<Type>> newBulkScheme
        (
            const fvMesh& mesh,
            const surfaceScalarField* faceFluxPtr,
            Istream& is
        );

        //- No copy construct
        momentumConsistentKuhn(const momentumConsistentKuhn&) = delete;

        //- No copy assignment
        void operator=(const momentumConsistentKuhn&) = delete;


public:

    //- Runtime type information
    TypeName("momentumConsistentKuhn");


    // Constructors

        //- Construct from mesh, faceFlux and Istream
        momentumConsistentKuhn
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
        momentumConsistentKuhn
        (
            const fvMesh& mesh,
            Istream& schemeData
        )
        :
            upwind<Type>(mesh, schemeData),
            tbulkScheme_(newBulkScheme(mesh, nullptr, schemeData))
        {}


    // Member Functions

        //- Blended interpolation weights (implicit part)
        virtual tmp<surfaceScalarField> weights
        (
            const GeometricField<Type, fvPatchField, volMesh>& vf
        ) const;

        //- Always corrected: the band carries the Kuhn correction
        virtual bool corrected() const
        {
            return true;
        }

        //- Band correction towards the Kuhn Favre value; bulk-scheme
        //  correction outside the band
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
