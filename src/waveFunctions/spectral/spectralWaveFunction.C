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

#include "spectralWaveFunction.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
  defineTypeNameAndDebug(spectralWaveFunction, 0);
  addToRunTimeSelectionTable
  (
   waveFunction,
   spectralWaveFunction,
   dict
   );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

/*Foam::implicitFunctions::spectralWaveFunction::spectralWaveFunction
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


Foam::spectralWaveFunction::spectralWaveFunction
(
     const dictionary& dict
)
  :
  amplitude_(dict.getOrDefault<scalar>("amplitude", 1)),
  waterLevel_(dict.getOrDefault<scalar>("waterLevel", 0)),
  phases_(dict.get<List<scalar>>("phases")),
  comps_(dict.get<List<scalar>>("comps")),
  kx_(dict.get<List<scalar>>("kx")),
  ky_(dict.get<List<scalar>>("ky"))
{
  //Info<<"comps="<<comps_<<endl;
}

// ************************************************************************* //
