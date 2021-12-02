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
    Foam::implicitFunctions::JONSWAPWaveFunction

Description

SourceFiles
    JONSWAPWaveFunction.C

\*---------------------------------------------------------------------------*/

#ifndef waveFunction_JONSWAPWaveFunction_H
#define waveFunction_JONSWAPWaveFunction_H

#include "spectralWaveFunction.H"
#include "mathematicalConstants.H"
#include "symmTensor2D.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
  //namespace implicitFunctions
  //{

/*---------------------------------------------------------------------------*\
                     Class JONSWAPWaveFunction Declaration
\*---------------------------------------------------------------------------*/

class JONSWAPWaveFunction
:
    public spectralWaveFunction
{
private:

    // Private Data
  const scalar amplitude_;
  const scalar waterLevel_;
  const scalar Tp_;
  const scalar Hs_;
  const int Nkx_;
  const int Nky_;
  const scalar gamma_;
  const scalar beta_;
  const scalar alphaJ_;
  const scalar Lx_;
  const scalar Ly_;
  
  /*const List<scalar> phases_;
  const List<scalar> comps_;
  const List<scalar> kx_;
  const List<scalar> ky_;*/

public:

    //- Runtime type information
    TypeName("JONSWAPWaveFunction");


    // Constructors

        //- Construct from components
        /*JONSWAPWaveFunction
        (
	 const scalar amplitude,
	 const scalar waterLevel,
	 const List<scalar>& phases,
	 const List<scalar>& comps,
	 const List<scalar>& kx,
	 const List<scalar>& ky
	 );*/ 

        //- Construct from dictionary (used by implicitFunctions)
  JONSWAPWaveFunction(const dictionary& dict);


    //- Destructor
  virtual ~JONSWAPWaveFunction() = default;


  virtual word name() const
  {
    char buffer [50];
    sprintf (buffer, "Hs%.3f_Tp%.3f_beta%.3f_gamma%.3f_Nk%d",Hs_,Tp_,beta_,gamma_,Nkx_);
    ///    word wname="X
    return buffer;
  }

  
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
