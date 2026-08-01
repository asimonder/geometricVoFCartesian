/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019 OpenCFD Ltd.
    Copyright (C) 2019-2020 DLR
    Copyright (C) 2025 Asim Önder 
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

#include "solitaryWaveLaitoneImplicitFunction.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace implicitFunctions
    {
        defineTypeNameAndDebug(solitaryWaveLaitoneImplicitFunction, 0);
        addToRunTimeSelectionTable
        (
            implicitFunction,
            solitaryWaveLaitoneImplicitFunction,
            dict
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::implicitFunctions::solitaryWaveLaitoneImplicitFunction::solitaryWaveLaitoneImplicitFunction
(
    const scalar height,
    const scalar depth,
    const scalar x0,
    const vector& direction,
    const vector& up,
    const vector& origin
)
:
    height_(height),
    depth_(depth),
    x0_(x0),
    direction_(normalised(direction)),
    up_(normalised(up)),
    origin_(origin)
{}


Foam::implicitFunctions::solitaryWaveLaitoneImplicitFunction::solitaryWaveLaitoneImplicitFunction
(
    const dictionary& dict
)
:
    height_(dict.get<scalar>("height")),
    depth_(dict.get<scalar>("depth")),
    x0_(dict.get<scalar>("x0")),
    direction_(dict.get<vector>("direction").normalise()),
    up_(dict.get<vector>("up").normalise()),
    origin_(dict.get<vector>("origin"))
{}


// ************************************************************************* //
