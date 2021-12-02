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

#include "Stokes3rdOrderWaveFunction.H"
#include "addToRunTimeSelectionTable.H"
#include "Pstream.H"
#include "Random.H"
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
  defineTypeNameAndDebug(Stokes3rdOrderWaveFunction, 0);
  addToRunTimeSelectionTable
  (
   waveFunction,
   Stokes3rdOrderWaveFunction,
   dict
   );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
Foam::Stokes3rdOrderWaveFunction::Stokes3rdOrderWaveFunction
(
     const dictionary& dict
)
  :
  amplitude_(dict.getOrDefault<scalar>("amplitude", 1)),
  waterLevel_(dict.getOrDefault<scalar>("waterLevel",0)),
  lamb_(dict.get<scalar>("lamb")),
  ak_(dict.get<scalar>("ak"))
{
  const scalar pi=constant::mathematical::pi;
  const scalar g=9.81;
   
  const label Nwaves=4; 
  kx_.setSize(Nwaves,0.0);
  ky_.setSize(Nwaves,0.0);
  comps_.setSize(Nwaves,0.0);
  phases_.setSize(Nwaves,0.0);

  for(int i=1;i<Nwaves;i++)
    {
      kx_[i]=2.0*pi*i/lamb_; 
      phases_[i]=0; 
    }

  scalar a=ak_/kx_[1];
  comps_[1]=a*(1-1./16*Foam::sqr(ak_));
  comps_[2]=a*(0.5*ak_);
  comps_[3]=a*(0.375*sqr(ak_));

}

// ************************************************************************* //
