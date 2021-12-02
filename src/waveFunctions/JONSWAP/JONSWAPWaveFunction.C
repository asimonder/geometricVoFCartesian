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

#include "JONSWAPWaveFunction.H"
#include "addToRunTimeSelectionTable.H"
#include "Pstream.H"
#include "Random.H"
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
  defineTypeNameAndDebug(JONSWAPWaveFunction, 0);
  addToRunTimeSelectionTable
  (
   waveFunction,
   JONSWAPWaveFunction,
   dict
   );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
Foam::JONSWAPWaveFunction::JONSWAPWaveFunction
(
     const dictionary& dict
)
  :
  amplitude_(dict.getOrDefault<scalar>("amplitude", 1)),
  waterLevel_(dict.getOrDefault<scalar>("waterLevel",0)),
  Hs_(dict.get<scalar>("Hs")),
  Tp_(dict.get<scalar>("Tp")),
  Nkx_(dict.get<int>("Nkx")),
  Nky_(dict.getOrDefault<int>("Nky",Nkx_)),
  gamma_(dict.get<scalar>("gamma")),
  beta_(dict.get<scalar>("beta")),
  alphaJ_(dict.get<scalar>("alphaJ")),
  //phi_(dict.get<scalar>("phi")),
  Lx_(dict.get<scalar>("Lx")),
  Ly_(dict.getOrDefault<scalar>("Ly",Lx_))
{
  const scalar pi=constant::mathematical::pi;
  const scalar g=9.81;
  bool is2D=(beta_==0?true:false);  
  label seed=0;
  if (Pstream::master())
    seed=(unsigned) time(0);
  Foam::reduce(seed, maxOp<label>());
  Info<<"seed="<<(unsigned) time(0)<<endl;
  Random rndPhase(seed); 

  const scalar omega_p=2*pi/Tp_;

  if (!is2D)
    {
      const scalar dkx=2*pi/Lx_;
      const scalar dky=2*pi/Ly_;
      //const scalar T1=0.834*Tp_;
 
      const label Nwaves=(Nkx_+1)*(2*Nky_+1); 
      kx_.setSize(Nwaves,0.0);
      ky_.setSize(Nwaves,0.0);
      comps_.setSize(Nwaves,0.0);
      phases_.setSize(Nwaves,0.0);

      label iC=0;
      scalar sigma=0;
      for(int ikx=0;ikx<Nkx_+1;ikx++)
	{
	  for(int iky=-Nky_;iky<Nky_+1;iky++)
	    {
	      if (ikx==0 && iky==0)
		continue;
	      
	      scalar kx=ikx*dkx;
	      scalar ky=iky*dky;
	      scalar k=Foam::sqrt(Foam::sqr(kx)+Foam::sqr(ky));
	      scalar theta=Foam::acos(kx/k);
	      
	      //scalar F=2.0/pi*Foam::sqr(Foam::cos(theta));
	      //scalar h=Foam::pow(gamma_,Foam::exp(-0.5*Foam::sqr((0.191*omega*T1-1)/sigma)));
	      //scalar Gw=alphaJ_*Foam::sqr(Hs_)*Foam::pow(T1,-4.)*Foam::pow(omega,-5.)*Foam::exp(-944./Foam::pow(T1*omega,4.0))*h;
	      if (abs(theta)<=beta_)
		{
		  scalar omega=Foam::sqrt(g*k);
		  scalar sigma=0.07;
		  if (omega>omega_p)
		    sigma=0.09;
		  scalar h=Foam::pow(gamma_,Foam::exp(-1.0/2*Foam::sqr((omega-omega_p)/(sigma*omega_p))));
		  scalar Gw=alphaJ_*Foam::sqr(Hs_)*Foam::pow(omega_p,4.)*Foam::pow(omega,-5)*Foam::exp(-5/4.*Foam::pow(omega/omega_p,-4))*h;
		  kx_[iC]=kx;
		  ky_[iC]=ky;
		  if (beta_==0)
		    comps_[iC]=Foam::sqrt(2*Foam::sqrt(g/k)*Gw*dkx);
		  else
		    {
		      scalar F=1.0/beta_*Foam::sqr(Foam::cos(pi/(2.0*beta_)*theta));
		      comps_[iC]=Foam::sqrt(2*Foam::sqrt(g/Foam::pow(k,3.))*Gw*F*dkx*dky);
		    }
		  phases_[iC]=rndPhase.globalPosition<scalar>(0,2*pi); 
		}
	      iC+=1;
	    }
	}
    }
  else
    {
      const scalar dkx=2*pi/Lx_;
      const label Nwaves=(Nkx_); 
      kx_.setSize(Nwaves,0.0);
      ky_.setSize(Nwaves,0.0);
      comps_.setSize(Nwaves,0.0);
      phases_.setSize(Nwaves,0.0);
      label iC=0;
      scalar sigma=0;
      for(int ikx=1;ikx<Nkx_+1;ikx++)
	{
	      scalar kx=ikx*dkx;
	      scalar omega=Foam::sqrt(g*kx);
	      scalar sigma=0.07;
	      if (omega>omega_p)
		sigma=0.09;
	      scalar h=Foam::pow(gamma_,Foam::exp(-1.0/2*Foam::sqr((omega-omega_p)/(sigma*omega_p))));
	      scalar Gw=alphaJ_*Foam::sqr(Hs_)*Foam::pow(omega_p,4.)*Foam::pow(omega,-5)*Foam::exp(-5/4.*Foam::pow(omega/omega_p,-4))*h;
	      kx_[ikx-1]=kx;
	      comps_[ikx-1]=Foam::sqrt(2*Foam::sqrt(g/kx)*Gw*dkx);
	      phases_[ikx-1]=rndPhase.globalPosition<scalar>(0,2*pi); 
	}
    }
}

// ************************************************************************* //
