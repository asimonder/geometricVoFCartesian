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
    Foam::implicitFunctions::Stokes3rdOrderWaveFunction

Description

SourceFiles
    Stokes3rdOrderWaveFunction.C

\*---------------------------------------------------------------------------*/

#ifndef waveFunction_Stokes3rdOrderWaveFunction_H
#define waveFunction_Stokes3rdOrderWaveFunction_H

#include "spectralWaveFunction.H"
#include "mathematicalConstants.H"
#include "symmTensor2D.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
  //namespace implicitFunctions
  //{

/*---------------------------------------------------------------------------*\
                     Class Stokes3rdOrderWaveFunction Declaration
\*---------------------------------------------------------------------------*/

class Stokes3rdOrderWaveFunction
:
    public spectralWaveFunction
{
private:

    // Private Data
  const scalar amplitude_;
  const scalar waterLevel_;
  const scalar lamb_;
  const scalar ak_;
  
  /*const List<scalar> phases_;
  const List<scalar> comps_;
  const List<scalar> kx_;
  const List<scalar> ky_;*/

public:

    //- Runtime type information
    TypeName("Stokes3rdOrderWaveFunction");


    // Constructors

        //- Construct from components
        /*Stokes3rdOrderWaveFunction
        (
	 const scalar amplitude,
	 const scalar waterLevel,
	 const List<scalar>& phases,
	 const List<scalar>& comps,
	 const List<scalar>& kx,
	 const List<scalar>& ky
	 );*/ 

        //- Construct from dictionary (used by implicitFunctions)
  Stokes3rdOrderWaveFunction(const dictionary& dict);


    //- Destructor
  virtual ~Stokes3rdOrderWaveFunction() = default;


  virtual word name() const
  {
    char buffer [50];
    sprintf (buffer, "ak%.3f",ak_);
    return buffer;
  }

  
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
