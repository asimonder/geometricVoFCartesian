/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019 OpenCFD Ltd.
    Copyright (C) 2019-2020 DLR
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

#include "flowerImplicitFunction.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace implicitFunctions
    {
        defineTypeNameAndDebug(flowerImplicitFunction, 0);
        addToRunTimeSelectionTable
        (
            implicitFunction,
            flowerImplicitFunction,
            dict
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::implicitFunctions::flowerImplicitFunction::flowerImplicitFunction
(
 const point& origin,
 const scalar radius,
 const scalar a,
 const scalar b,
 const scalar p
)
:
    origin_(origin),
    radius_(radius),
    a_(a),
    b_(b),
    p_(p)
    //direction_(normalised(direction)),
    //project_(tensor::I - direction_*direction_) // outer product
{}


/*Foam::implicitFunctions::flowerImplicitFunction::flowerImplicitFunction
(
    const dictionary& dict
)
:
    flowerImplicitFunction
    (
        dict.get<point>("origin"),
        dict.get<scalar>("radius"),
        dict.getOrDefault<scalar>("scale", 1),
        dict.get<vector>("direction")
    )
    {}*/

Foam::implicitFunctions::flowerImplicitFunction::flowerImplicitFunction
 (
     const dictionary& dict
 )
 :
     // __INTEL_COMPILER bug with inheriting constructors?? (issue #1821)
     origin_(dict.get<point>("origin")),
     radius_(dict.get<scalar>("radius")),
     a_(dict.get<scalar>("a")),
     b_(dict.get<scalar>("b")),
     p_(dict.get<scalar>("p"))
     //scale_(dict.getOrDefault<scalar>("scale", 1)),
     //direction_(dict.get<vector>("direction").normalise()),
     //project_(tensor::I - direction_*direction_) // outer product
 {}
  
// ************************************************************************* //
