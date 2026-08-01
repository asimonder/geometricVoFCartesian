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
    Foam::isoAdvector

Description
    geometricAdvector wrapper around the stock isoAdvection class
    (Roenby, Bredmose & Jasak 2016), selected as

        advector isoAdvector;

    Source-free advection only; the consistent per-phase momentum
    interface is not provided by isoAdvector and fails fatally
    through the base class.

SourceFiles
    isoAdvector.C

\*---------------------------------------------------------------------------*/

#ifndef isoAdvectorGA_H
#define isoAdvectorGA_H

#include "geometricAdvector.H"
#include "isoAdvection.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class isoAdvector Declaration
\*---------------------------------------------------------------------------*/

class isoAdvector
:
    public geometricAdvector
{
    // Private data

        isoAdvection advector_;


public:

    //- Runtime type information
    TypeName("isoAdvector");


    // Constructors

        isoAdvector
        (
            volScalarField& alpha1,
            const surfaceScalarField& phi,
            const volVectorField& U
        )
        :
            advector_(alpha1, phi, U)
        {}


    //- Destructor
    virtual ~isoAdvector() = default;


    // Member Functions

        virtual void advect()
        {
            zeroField Su;
            zeroField Sp;
            advector_.advect(Sp, Su);
        }

        virtual const surfaceScalarField& alphaPhi() const
        {
            return advector_.alphaPhi();
        }

        virtual reconstructionSchemes& surf()
        {
            return advector_.surf();
        }
};


} // End namespace Foam

#endif

// ************************************************************************* //
