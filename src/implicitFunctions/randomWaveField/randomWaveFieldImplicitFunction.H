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
    Foam::implicitFunctions::randomWaveFieldImplicitFunction

Description

SourceFiles
    randomWaveFieldImplicitFunction.C

\*---------------------------------------------------------------------------*/

#ifndef implicitFunction_randomWaveFieldImplicitFunction_H
#define implicitFunction_randomWaveFieldImplicitFunction_H

#include "implicitFunction.H"
#include "mathematicalConstants.H"
#include "symmTensor2D.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace implicitFunctions
{

/*---------------------------------------------------------------------------*\
                     Class randomWaveFieldImplicitFunction Declaration
\*---------------------------------------------------------------------------*/

class randomWaveFieldImplicitFunction
:
    public implicitFunction
{
private:

    // Private Data
  const scalar amplitude_;
  const scalar waterLevel_;

protected:
  List<scalar> phases_;
  List<scalar> comps_;
  List<scalar> kx_;
  List<scalar> ky_;

public:

  //- Runtime type information
  TypeName("randomWaveField");

    // Constructors
  /*randomWaveFieldImplicitFunction()
    :
    amplitude_(1.0),
    waterLevel_(0)
  {
    //amplitude_=1;
    // waterLevel_=0;
    };*/
  
  //- Construct from components
  randomWaveFieldImplicitFunction
  (
   const scalar amplitude,
   const scalar waterLevel,
   const List<scalar>& phases,
   const List<scalar>& comps,
   const List<scalar>& kx,
   const List<scalar>& ky
   );
  
  //- Construct from dictionary (used by implicitFunctions)
  explicit randomWaveFieldImplicitFunction(const dictionary& dict);


    //- Destructor
  virtual ~randomWaveFieldImplicitFunction() = default;


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
    
    return H-z-waterLevel_;
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

    G.z()=-1;
    
    return G;
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
    return H;
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

  void setWaterLevel(scalar H)
  {
    waterLevel_=H;
  }
  
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace implicitFunctions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
