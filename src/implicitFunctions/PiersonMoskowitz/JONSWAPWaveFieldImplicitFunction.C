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
  XFetch_(dict.get<scalar>("XFetch")),
  UWind_(dict.get<scalar>("UWind")),
  omegaMin_(dict.lookupOrDefault<scalar>("omegaMin",0.5)),
  omegaMax_(dict.lookupOrDefault<scalar>("omegaMax",3)),
  NOmega_(dict.lookupOrDefault<label>("NOmega",26)),
  NTheta_(dict.lookupOrDefault<label>("NTheta",37)), 
  gamma_(dict.lookupOrDefault<scalar>("gamma",3.3))
  //  phases_(dict.get<List<scalar>>("phases")),
  //comps_(dict.get<List<scalar>>("comps")),
  //kx_(dict.get<List<scalar>>("kx")),
  //ky_(dict.get<List<scalar>>("ky"))
{
  //Info<<"comps="<<comps_<<endl;

  const scalar pi=constant::mathematical::pi;
  
  label seed=0;
  if (Pstream::master())
    seed=(unsigned) time(0);
  Foam::reduce(seed, maxOp<label>());
  Info<<"seed="<<(unsigned) time(0)<<endl;
  Random rndPhase(seed); 


  const scalar g=9.81;
  //const scalar gamma=3.3;
  const scalar ksi=XFetch_*g/Foam::sqr(UWind_);
  const scalar alpha=0.076*Foam::pow(ksi,-0.22);
  const scalar omega_p=7.0*pi*(g/UWind_)*Foam::pow(ksi,-0.33);

  List<scalar> omegaL(NOmega_,0.0);
  List<scalar> thetaL(NTheta_,0.0);

  const scalar dTheta=pi/(NTheta_-1);
  const scalar dOmega=omega_p*(omegaMax_-omegaMin_)/(NOmega_-1);

  label iC=0;
  forAll(omegaL,iOmega)
    omegaL[iOmega]=omegaMin_+iOmega*dOmega;

  forAll(thetaL,iTheta)
    thetaL[iTheta]=-pi/2+iTheta*dTheta;

  
  label Nwaves=NTheta_*NOmega_; 
  kx_.setSize(Nwaves,0.0);
  ky_.setSize(Nwaves,0.0);
  comps_.setSize(Nwaves,0.0);
  phases_.setSize(Nwaves,0.0);

  iC=0;
  scalar sigma=0;
  forAll(omegaL,iOmega)
    {
      if (omegaL[iOmega]<omega_p)
	sigma=0.07;
      else
	sigma=0.09;
      const scalar h=Foam::pow(gamma_,Foam::exp(-0.5*Foam::sqr((omegaL[iOmega]/omega_p-1)/sigma)));
      const scalar G=alpha*Foam::sqr(g)*Foam::pow(omegaL[iOmega],-5.)*Foam::exp(-5/4.*Foam::pow(omega_p/omegaL[iOmega],4.0))*h;
      const scalar k=Foam::sqr(omegaL[iOmega])/g;
      forAll(thetaL,iTheta)
	{
	  const scalar F=2/pi*Foam::sqr(Foam::cos(thetaL[iTheta]));
	  kx_[iC]=k*Foam::cos(thetaL[iTheta]);
	  ky_[iC]=k*Foam::sin(thetaL[iTheta]);
	  comps_[iC]=Foam::sqrt(2.*G*F*dTheta*dOmega);
	  phases_[iC]=rndPhase.globalPosition<scalar>(0,2*pi); 
	  iC+=1;
	}
    }
}

// ************************************************************************* //
