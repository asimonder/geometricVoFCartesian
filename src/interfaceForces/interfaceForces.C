/*---------------------------------------------------------------------------* \
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "interfaceForces.H"
//#include "alphaContactAngleTwoPhaseFvPatchScalarField.H"
#include "mathematicalConstants.H"
#include "surfaceInterpolate.H"
#include "fvcDiv.H"
//#include "fvcReconstruct.H"
#include "fvcGrad.H"
#include "fvcSnGrad.H"
#include "unitConversion.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Correction for the boundary condition on the unit normal nHat on
// walls to produce the correct contact angle.

// The dynamic contact angle is calculated from the component of the
// velocity on the direction of the interface, parallel to the wall.



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::interfaceForces::interfaceForces
(
    const volScalarField& alpha1,
    const volVectorField& U,
    const IOdictionary& dict
)
:
    transportPropertiesDict_(dict),

    sigmaPtr_(surfaceTensionModel::New(dict,alpha1.mesh())),

    curvature_(curvatureModel::New(dict,alpha1,U)),

    alpha1_(alpha1),
    U_(U)
{
  //calculateK();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

/*Foam::tmp<Foam::volScalarField>
Foam::interfaceForces::sigmaK() const
{
    return sigmaPtr_->sigma()*K_;
    }*/

Foam::tmp<Foam::volScalarField>
Foam::interfaceForces::sigma() const
{
    return sigmaPtr_->sigma();
}

Foam::tmp<Foam::volScalarField>
Foam::interfaceForces::K() const
{
    return curvature_->K();
}

Foam::tmp<Foam::surfaceScalarField>
Foam::interfaceForces::Kf() const
{
    return curvature_->Kf();
}

Foam::tmp<Foam::surfaceScalarField>
Foam::interfaceForces::surfaceTensionForce() const
{
  //Foam::surfaceScalarField Kf(fvc::interpolate(K()*0.0));
  /*Foam::surfaceScalarField Kf(fvc::interpolate(K()));
  const Foam::volScalarField& Kc=K();
  const Foam::labelList& nei=alpha1_.mesh().faceNeighbour();
  const Foam::labelList& own=alpha1_.mesh().faceOwner();
  forAll(Kf,iFace)
    {
      scalar K1=Kc[own[iFace]];      
      scalar K2=0;
      if (nei[iFace]>-1)
	{
	  K2=Kc[nei[iFace]];
	}
      else
	{
	  K2=-1.;
	  Info<<"face "<<iFace<<" is on the boundary..."<<endl;
	}
      if (K1*K2==0 and Kf[iFace]!=0) 
	Kf[iFace]*=2;
    }
    return fvc::interpolate(sigma())*Kf*fvc::snGrad(alpha1_);*/
  return fvc::interpolate(sigma())*Kf()*fvc::snGrad(alpha1_);
  // return fvc::interpolate(sigma()*K())*fvc::snGrad(alpha1_);
  //return fvc::snGrad(sigma()*K()*alpha1_);
}

/*Foam::tmp<Foam::volVectorField>
Foam::interfaceForces::surfaceTensionForce() const
{
  //return  sigma()*K()*fvc::reconstruct(fvc::snGrad(alpha1_));
  volVectorField sigmaForce(K()*sigma()*fvc::reconstruct(fvc::snGrad(alpha1_)));
  //forAll(sigmaForce,iCell)
  // {
  //   sigmaForce[iCell]*=K().values()[iCell];//*sigma();
  // }
  return  sigmaForce; //fvc::reconstruct(fvc::snGrad(alpha1_));
  }*/

Foam::tmp<Foam::volScalarField>
Foam::interfaceForces::nearInterface() const
{
    return pos0(alpha1_ - 0.01)*pos0(0.99 - alpha1_);
}


bool Foam::interfaceForces::read()
{
  //alpha1_.mesh().solverDict(alpha1_.name()).readEntry("cAlpha", cAlpha_);
    sigmaPtr_->readDict(transportPropertiesDict_);

    return true;
}

void Foam::interfaceForces::correct()
{
     return curvature_->correct();
}
// ************************************************************************* //
