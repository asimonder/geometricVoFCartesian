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

#include "JONSWAPWaveFieldImplicitFunction.H"
#include "addToRunTimeSelectionTable.H"
#include "Pstream.H"
#include "Random.H"
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace implicitFunctions
    {
        defineTypeNameAndDebug(JONSWAPWaveFieldImplicitFunction, 0);
        addToRunTimeSelectionTable
        (
            implicitFunction,
            JONSWAPWaveFieldImplicitFunction,
            dict
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

/*Foam::implicitFunctions::JONSWAPWaveFieldImplicitFunction::JONSWAPWaveFieldImplicitFunction
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


Foam::implicitFunctions::JONSWAPWaveFieldImplicitFunction::JONSWAPWaveFieldImplicitFunction
(
     const dictionary& dict
)
  :
  amplitude_(dict.getOrDefault<scalar>("amplitude", 1)),
  waterLevel_(dict.getOrDefault<scalar>("waterLevel",0)),
  Tp_(dict.get<scalar>("Hs")),
  Hs_(dict.get<scalar>("Tp")),
  Nk_(dict.get<scalar>("Nk")),
  gamma_(dict.get<scalar>("gamma")),
  beta_(dict.get<scalar>("beta")),
  alpha_(dict.get<scalar>("alpha")),
  phi_(dict.get<scalar>("phi")),
  L_(dict.getOrDefault<scalar>("domainLength",1)),
{
  const scalar pi=constant::mathematical::pi;
  const scalar g=9.81;
  
  label seed=0;
  if (Pstream::master())
    seed=(unsigned) time(0);
  Foam::reduce(seed, maxOp<label>());
  Info<<"seed="<<(unsigned) time(0)<<endl;
  Random rndPhase(seed); 

  const scalar dk=2*pi/L_;
  const scalar T1=0.834*Tp;
  const scalar omega_p=2*pi/Tp
  
 
  label Nwaves=(Nk_+1)*(Nk_+1); 
  kx_.setSize(Nwaves,0.0);
  ky_.setSize(Nwaves,0.0);
  comps_.setSize(Nwaves,0.0);
  phases_.setSize(Nwaves,0.0);

  iC=0;
  scalar sigma=0;
  for(int ikx=0;ikx<Nkx+1;ikx++)
    {
      for(int iky=-Nky;iky<Nky+1;iky++)
	{
	  if (ikx==0 && iky==0)
            continue;

	  scalar kx=ikx*dkx;
	  scalar ky=iky*dky;
	  scalar k=Foam::sqrt(Foam:sqr(kx)+Foam::sqr(ky));
	  scalar theta=Foam::acos(kx/k);
	  scalar omega=Foam::sqrt(g*k);
	  scalar sigma=0.07;
	  if (omega>omega_p)
	    sigma=0.09;

	  scalar F=2.0/pi*Foam::sqr(Foam::cos(theta));
	  scalar h=Foam::pow(gamma_,Foam::exp(-0.5*Foam::sqr((0.191*omega*T1-1)/sigma)));
	  scalar Gw=alpha_*Foam::sqr(Hs)*Foam::pow(T1,-4.)*Foam::pow(omega,-5.)*Foam::exp(-944./Foam::pow(T1*omega,4.0))*h;
	  if (Foam::abs(theta)<=beta)
	    {
	      kx_[iC]=kx
	      ky_[iC]=ky;
	      comps_[iC]=Foam::sqrt(Foam::sqrt(g/Foam::pow(k,3.)*Gw*F*dkx*dky);
	      phases_[iC]=rndPhase.globalPosition<scalar>(0,2*pi); 
	    }

	  iC+=1;
	}
    }

}

// ************************************************************************* //
