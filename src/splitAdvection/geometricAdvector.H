/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2026 Asim Onder
     \\/     M anipulation  |
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
    Foam::geometricAdvector

Description
    Runtime-selectable interface for geometric VOF advection so that
    solvers are written once against this surface and the advector is
    chosen per case from fvSolution:

        "alpha.*"
        {
            advector        wy;          // or isoAdvector (default)
            ...
        }

    The consistent per-phase momentum transport advect(q1, q2) is
    optional for an implementation; the base version fails fatally so
    a case cannot silently run momentum-inconsistent transport where
    consistency was requested.

SourceFiles
    geometricAdvector.C

\*---------------------------------------------------------------------------*/

#ifndef geometricAdvector_H
#define geometricAdvector_H

#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "reconstructionSchemes.H"
#include "runTimeSelectionTables.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class geometricAdvector Declaration
\*---------------------------------------------------------------------------*/

class geometricAdvector
{
public:

    //- Runtime type information
    TypeName("geometricAdvector");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        geometricAdvector,
        components,
        (
            volScalarField& alpha1,
            const surfaceScalarField& phi,
            const volVectorField& U
        ),
        (alpha1, phi, U)
    );


    // Constructors

        geometricAdvector() = default;

        //- Select from the "advector" keyword of the alpha solver
        //  dictionary (default isoAdvector)
        static autoPtr<geometricAdvector> New
        (
            volScalarField& alpha1,
            const surfaceScalarField& phi,
            const volVectorField& U
        );


    //- Destructor
    virtual ~geometricAdvector() = default;


    // Member Functions

        //- Advance alpha1 by one time step
        virtual void advect() = 0;

        //- Advance alpha1 and the phase momenta consistently
        virtual void advect(volVectorField& q1, volVectorField& q2)
        {
            FatalErrorInFunction
                << type() << " does not implement consistent "
                << "per-phase momentum transport"
                << exit(FatalError);
        }

        //- Liquid volume face flux of the completed step [m3/s]
        virtual const surfaceScalarField& alphaPhi() const = 0;

        //- Access the interface reconstruction
        virtual reconstructionSchemes& surf() = 0;

        //- Step diagnostics (implementations may override)
        virtual label nClamped() const { return 0; }
        virtual scalar maxExcursion() const { return 0; }
};


} // End namespace Foam

#endif

// ************************************************************************* //
