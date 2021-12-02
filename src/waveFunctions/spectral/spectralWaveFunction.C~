/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 Asim Onder
-------------------------------------------------------------------------------
License
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

#include "randomWaveFieldImplicitFunction.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace implicitFunctions
    {
        defineTypeNameAndDebug(randomWaveFieldImplicitFunction, 0);
        addToRunTimeSelectionTable
        (
            implicitFunction,
            randomWaveFieldImplicitFunction,
            dict
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

/*Foam::implicitFunctions::randomWaveFieldImplicitFunction::randomWaveFieldImplicitFunction
(
 const scalar amplitude,
 const scalar waterLevel,
 const List<scalar>& phases,
 const List<scalar>& comps,
 const List<scalar>& kx,
 const List<scalar>& ky
)
:
  amplitude_(amplitude),
  waterLevel_(waterLevel),
  phases_(phases),
  comps_(comps),
  kx_(kx),
  ky_(ky)
  {}*/


Foam::implicitFunctions::randomWaveFieldImplicitFunction::randomWaveFieldImplicitFunction
(
     const dictionary& dict
)
  :
  amplitude_(dict.getOrDefault<scalar>("amplitude", 1)),
  waterLevel_(dict.getOrDefault<scalar>("waterLevel", 0)),
  phases_(dict.getOrDefault<List<scalar>>("phases",0)),
  comps_(dict.getOrDefault<List<scalar>>("comps",0)),
  kx_(dict.getOrDefault<List<scalar>>("kx",0)),
  ky_(dict.getOrDefault<List<scalar>>("ky",0))
{
  //Info<<"comps="<<comps_<<endl;
}

// ************************************************************************* //
