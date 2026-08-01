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
    Foam::momentumConsistentFavreMinmod

Description
    Interface-band blended convection interpolation following Kuhn,
    Deskos, Sprague (2023), Computers and Fluids 252:105770, Eq. 35.
    Band faces (interfaceBand mesh object) are implicit pure upwind
    plus an explicit deferred correction towards the Favre-averaged
    minmod-limited face value

        U_f = interp_w(rho U) / interp_w(rho),

    with ONE shared weight field w from the scalar Minmod limiter of
    rho (stock vector limiters evaluate their ratio on magSqr, so
    per-field limiters break uniform-flow cancellation). The face
    value is convex in the two cell values for any w in [0,1], so the
    flux cannot remove more momentum than the cells hold (Kuhn Fig. 3
    argument). Uniform flow cancels the correction exactly. Bulk
    faces use an arbitrary scheme read from the fvSchemes line. Ramp
    rings of the band blend the corrections smoothly.

    A symmetric weights-only form (one shared weight, density-folded)
    was tested first and is unstable: with the light phase upwind the
    folded weight turns downwind. Hence the deferred-correction form.

    Deviations from Kuhn et al.: no half-step temporal extrapolation
    (their Eq. 31 second term; our flux is the interval average), no
    PPM transverse gradients (their Eq. 34), and the shared limiter
    weight described above instead of per-field limiters.

    Usage:
    \verbatim
    divSchemes
    {
        div(rhoPhi,U)  Gauss momentumConsistentFavreMinmod linear;
    }
    \endverbatim

    With no bulk scheme given, linear is used. A leading numeric token
    (legacy interfaceThreshold argument) is accepted and ignored.

\*---------------------------------------------------------------------------*/

#include "upwind.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifndef momentumConsistentFavreMinmod_H
#define momentumConsistentFavreMinmod_H

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class momentumConsistentFavreMinmod Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class momentumConsistentFavreMinmod
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
        momentumConsistentFavreMinmod
        (
            const momentumConsistentFavreMinmod&
        ) = delete;

        //- No copy assignment
        void operator=(const momentumConsistentFavreMinmod&) = delete;


public:

    //- Runtime type information
    TypeName("momentumConsistentFavreMinmod");


    // Constructors

        //- Construct from mesh, faceFlux and Istream
        momentumConsistentFavreMinmod
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
        momentumConsistentFavreMinmod
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

        //- Always corrected: the band carries the Favre correction
        virtual bool corrected() const
        {
            return true;
        }

        //- Bulk-scheme correction, zeroed on interface-band faces
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
