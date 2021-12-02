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
    Foam::implicitFunctions::JONSWAPWaveFieldImplicitFunction

Description

SourceFiles
    JONSWAPWaveFieldImplicitFunction.C

\*---------------------------------------------------------------------------*/

#ifndef implicitFunction_JONSWAPWaveFieldImplicitFunction_H
#define implicitFunction_JONSWAPWaveFieldImplicitFunction_H

#include "randomWaveFieldImplicitFunction.H"
#include "mathematicalConstants.H"
#include "symmTensor2D.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace implicitFunctions
{

/*---------------------------------------------------------------------------*\
                     Class JONSWAPWaveFieldImplicitFunction Declaration
\*---------------------------------------------------------------------------*/

class JONSWAPWaveFieldImplicitFunction
:
    public randomWaveFieldImplicitFunction
{
private:

  // Private Data
  const scalar amplitude_;
  const scalar waterLevel_;
  const scalar Tp_;
  const scalar Hs_;
  const label Nk_;
  const scalar gamma_;
  const scalar beta_;
  const scalar alpha_;
  const scalar L_;
  
  /*const List<scalar> phases_;
  const List<scalar> comps_;
  const List<scalar> kx_;
  const List<scalar> ky_;*/

public:

    //- Runtime type information
    TypeName("JONSWAPWaveField");


    // Constructors

        //- Construct from components
        /*JONSWAPWaveFieldImplicitFunction
        (
	 const scalar amplitude,
	 const scalar waterLevel,
	 const List<scalar>& phases,
	 const List<scalar>& comps,
	 const List<scalar>& kx,
	 const List<scalar>& ky
	 );*/ 

        //- Construct from dictionary (used by implicitFunctions)
        JONSWAPWaveFieldImplicitFunction(const dictionary& dict);


    //- Destructor
  virtual ~JONSWAPWaveFieldImplicitFunction() = default;

  
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace implicitFunctions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
