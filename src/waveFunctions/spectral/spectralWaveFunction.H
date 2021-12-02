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

Class
    Foam::implicitFunctions::spectralWaveFunction

Description

SourceFiles
    spectralWaveFunction.C

\*---------------------------------------------------------------------------*/

#ifndef implicitFunction_spectralWaveFunction_H
#define implicitFunction_spectralWaveFunction_H

#include "waveFunction.H"
#include "mathematicalConstants.H"
//#include "symmTensor2D.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class spectralWaveFunction Declaration
\*---------------------------------------------------------------------------*/

class spectralWaveFunction
:
    public waveFunction
{
private:

protected:
  scalar amplitude_;
  scalar waterLevel_;

  List<scalar> phases_;
  List<scalar> comps_;
  List<scalar> kx_;
  List<scalar> ky_;

public:

  //- Runtime type information
  TypeName("spectralWaveFunction");

    // Constructors
  spectralWaveFunction()
    :
    amplitude_(1.0),
    waterLevel_(0.0)
  {
    phases_.setSize(1,0);
    comps_.setSize(1,0);
    kx_.setSize(1,0);
    ky_.setSize(1,0);
  };
  
  /*spectralWaveFunction()
    :
    amplitude_(1.0),
    waterLevel_(0)
  {
    //amplitude_=1;
    // waterLevel_=0;
    };*/
  
  //- Construct from components
  /*spectralWaveFunction
  (
   const scalar amplitude,
   const scalar waterLevel,
   const List<scalar>& phases,
   const List<scalar>& comps,
   const List<scalar>& kx,
   const List<scalar>& ky
   );*/
  
  //- Construct from dictionary (used by implicitFunctions)
  spectralWaveFunction(const dictionary& dict);


    //- Destructor
  virtual ~spectralWaveFunction() = default;


    // Member Functions

  virtual scalar value(const vector& p) const
  {
    const scalar x=p.x();
    const scalar y=p.y();
    const scalar z=p.z();
    
    scalar H=0;
    forAll (comps_,i)
      {
	//Info<<"comps[i]="<<comps_[i]<<endl;
	H+=comps_[i]*Foam::cos(kx_[i]*x+ky_[i]*y+phases_[i]);
      }
    
    return amplitude_*H-z-waterLevel_;
  }

  virtual vector grad(const vector& p) const
  {
    const scalar x=p.x();
    const scalar y=p.y();

    vector G(0,0,0);
    forAll (comps_,i)
      {
	G.x()-=comps_[i]*kx_[i]*Foam::sin(kx_[i]*x+ky_[i]*y+phases_[i]);
	G.y()-=comps_[i]*ky_[i]*Foam::sin(kx_[i]*x+ky_[i]*y+phases_[i]);
      }

    G.z()=-1.0/amplitude_;
    
    return amplitude_*G;
  }

  virtual scalar distanceToSurfaces(const vector& p) const
  {
    NotImplemented;
    return 0;
  }

  virtual symmTensor2D hessian(const vector& p) const
  {
    //NotImplemented;
    
    const scalar x=p.x();
    const scalar y=p.y();
    
    symmTensor2D H(0.0,0.0,0.0);
    forAll (comps_,i)
      {
	H.xx()+=-comps_[i]*kx_[i]*kx_[i]*Foam::cos(kx_[i]*x+ky_[i]*y+phases_[i]);
	H.yy()+=-comps_[i]*ky_[i]*ky_[i]*Foam::cos(kx_[i]*x+ky_[i]*y+phases_[i]);
	H.xy()+=-comps_[i]*kx_[i]*ky_[i]*Foam::cos(kx_[i]*x+ky_[i]*y+phases_[i]);
      }

    //H.yx()=H.xy();	  
    return amplitude_*H;
  }
  
  virtual scalar curvature(const vector& p) const
  {
    //const scalar x=p.x();
    //const scalar y=p.y();
    
    scalar Hx=grad(p).x();
    scalar Hy=grad(p).y();
    scalar Hxx=hessian(p).xx();
    scalar Hyy=hessian(p).yy();
    scalar Hxy=hessian(p).xy();
    
    /*    forAll (comps_,i)
      {
	Hx+=-comps_[i]*kx_[i]*Foam::sin(kx_[i]*x+ky_[i]*y+phases_[i]);
	Hy+=-comps_[i]*ky_[i]*Foam::sin(kx_[i]*x+ky_[i]*y+phases_[i]);
	Hxx+=-comps_[i]*kx_[i]*kx_[i]*Foam::cos(kx_[i]*x+ky_[i]*y+phases_[i]);
	Hyy+=-comps_[i]*ky_[i]*ky_[i]*Foam::cos(kx_[i]*x+ky_[i]*y+phases_[i]);
	Hxy+=-comps_[i]*kx_[i]*ky_[i]*Foam::cos(kx_[i]*x+ky_[i]*y+phases_[i]);
	}*/

    scalar kappa=(Hxx+Hyy+Hxx*Hy*Hy+Hyy*Hx*Hx-2.0*Hxy*Hx*Hy)/Foam::pow(1.0+Hx*Hx+Hy*Hy,1.5);
    
    return -kappa;
  }

  virtual word name() const
  {
    return "Spectral_Wave";
  }

  virtual scalar getWaterLevel()
  {
    return waterLevel_;
  }

  virtual void setWaterLevel(scalar H)
  {
    waterLevel_=H;
  }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
