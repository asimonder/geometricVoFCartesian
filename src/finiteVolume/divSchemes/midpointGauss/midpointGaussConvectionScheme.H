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

Class
    Foam::fv::midpointGaussConvectionScheme

Description
    Gauss convection scheme assembling scheme (S) of the note
    "Consistent second-order time integration for isoAdvector".

    The interval-averaged geometric mass flux (rhoPhi) multiplies the
    trapezoidal average of the transported field. One half is assembled
    implicitly at n+1 and one half explicitly from oldTime(), both halves
    sharing the interpolation weights evaluated once at the current level.
    The weights sum to one at both levels, so uniform translation is
    preserved exactly for any time step and density ratio.

    Usage (with stock Euler ddt):
        div(rhoPhi,U)  midpointGauss momentumConsistentLinear linear;

SourceFiles
    midpointGaussConvectionScheme.C

\*---------------------------------------------------------------------------*/

#ifndef midpointGaussConvectionScheme_H
#define midpointGaussConvectionScheme_H

#include "convectionScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fv
{

/*---------------------------------------------------------------------------*\
              Class midpointGaussConvectionScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class midpointGaussConvectionScheme
:
    public fv::convectionScheme<Type>
{
    // Private data

        tmp<surfaceInterpolationScheme<Type>> tinterpScheme_;


    // Private Member Functions

        //- No copy construct
        midpointGaussConvectionScheme
        (
            const midpointGaussConvectionScheme&
        ) = delete;

        //- No copy assignment
        void operator=(const midpointGaussConvectionScheme&) = delete;


public:

    //- Runtime type information
    TypeName("midpointGauss");


    // Constructors

        //- Construct from flux and interpolation scheme
        midpointGaussConvectionScheme
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            const tmp<surfaceInterpolationScheme<Type>>& scheme
        )
        :
            convectionScheme<Type>(mesh, faceFlux),
            tinterpScheme_(scheme)
        {}

        //- Construct from flux and Istream
        midpointGaussConvectionScheme
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            Istream& is
        )
        :
            convectionScheme<Type>(mesh, faceFlux),
            tinterpScheme_
            (
                surfaceInterpolationScheme<Type>::New(mesh, faceFlux, is)
            )
        {}


    // Member Functions

        const surfaceInterpolationScheme<Type>& interpScheme() const;

        //- Trapezoidal face value, 0.5*(interp(vf) + interp(vf.oldTime()))
        tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> interpolate
        (
            const surfaceScalarField&,
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const;

        tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> flux
        (
            const surfaceScalarField&,
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const;

        //- Half implicit at n+1, half explicit from oldTime(),
        //  shared current-level weights
        tmp<fvMatrix<Type>> fvmDiv
        (
            const surfaceScalarField&,
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const;

        tmp<GeometricField<Type, fvPatchField, volMesh>> fvcDiv
        (
            const surfaceScalarField&,
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "midpointGaussConvectionScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
