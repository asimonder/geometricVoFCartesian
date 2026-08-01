/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019 OpenCFD Ltd.
    Copyright (C) 2019-2020 DLR
    Copyright (C) 2025 Asim Önder 
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

Class
    Foam::implicitFunctions::solitaryWaveLaitoneImplicitFunction

Description

    Laitone’s second-order theory (Laitone, 1960)

SourceFiles
    solitaryWaveLaitoneImplicitFunction.C

\*---------------------------------------------------------------------------*/

#ifndef implicitFunction_solitaryWaveLaitoneImplicitFunction_H
#define implicitFunction_solitaryWaveLaitoneImplicitFunction_H

#include "implicitFunction.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace implicitFunctions
{

/*---------------------------------------------------------------------------*\
                     Class solitaryWaveLaitoneImplicitFunction Declaration
\*---------------------------------------------------------------------------*/

class solitaryWaveLaitoneImplicitFunction
:
    public implicitFunction
{
private:

    // Private Data

  //- Origin point
  //const scalar period_;

  //- Radius
  //const scalar phase_;
  
  //const scalar amplitude_;
  
  const vector up_;

  const vector direction_;

  const vector origin_;

      // Parameters read from dictionary
  const scalar height_;    // H
  const scalar depth_;     // d
  const scalar x0_;    // x0

public:

    //- Runtime type information
    TypeName("solitaryWaveLaitone");


    // Constructors

        //- Construct from components
  solitaryWaveLaitoneImplicitFunction
  (
   const scalar height,
   const scalar depth,
   const scalar x0,
   const vector& direction,
   const vector& up,
   const vector& origin
   );

        //- Construct from dictionary (used by implicitFunctions)
  explicit solitaryWaveLaitoneImplicitFunction(const dictionary& dict);


    //- Destructor
  virtual ~solitaryWaveLaitoneImplicitFunction() = default;

  // Member Functions
  virtual scalar value(const vector& p) const
  {
    const scalar x = (p - origin_) & direction_;
    const scalar z = (p - origin_) & up_;

    const scalar xLoc = sqrt(3.0 * height_ / (4.0 * depth_))
                   * (1.0 - 5.0 * height_ / (8.0 * depth_))
                   * (x - x0_) / depth_;

    const scalar sechX = 1.0 / cosh(xLoc);
    const scalar sech2X = sechX * sechX;

    const scalar eta =
        height_ * sech2X
      - 0.75 * sqr(height_ / depth_) * depth_ * sech2X * (1.0 - sech2X);

    return -(z - eta);
  }
  
  virtual vector grad(const vector& p) const
  {
    NotImplemented;
    return Zero;
  }

  virtual scalar distanceToSurfaces(const vector& p) const
  {
    NotImplemented;
    return 0;
  }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace implicitFunctions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
