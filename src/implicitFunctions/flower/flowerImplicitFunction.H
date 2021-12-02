/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019 OpenCFD Ltd.
    Copyright (C) 2019-2020 DLR
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
    Foam::implicitFunctions::flowerImplicitFunction

Description
    creates a infintite long flowerImplicitFunction

    Original code supplied by Henning Scheufler, DLR (2019)

SourceFiles
    flowerImplicitFunction.C

\*---------------------------------------------------------------------------*/

#ifndef implicitFunction_flowerImplicitFunction_H
#define implicitFunction_flowerImplicitFunction_H

#include "implicitFunction.H"
#include "point.H"
#include "tensor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace implicitFunctions
{

/*---------------------------------------------------------------------------*\
                   Class flowerImplicitFunction Declaration
\*---------------------------------------------------------------------------*/

class flowerImplicitFunction
:
    public implicitFunction
{
private:

    // Private Data

        //- Origin point
        const point origin_;

        //- Radius
        const scalar radius_;

  const scalar a_;
  const scalar b_;
  const scalar p_;
  //const vector direction_;

  //const tensor project_;


public:

    //- Runtime type information
    TypeName("flower");


    // Constructors

        //- Construct from components
        flowerImplicitFunction
        (
            const point& origin,
            const scalar radius,
            const scalar a,
	    const scalar b,
	    const scalar p
        );

        //- Construct from dictionary
        explicit flowerImplicitFunction(const dictionary& dict);


    //- Destructor
    virtual ~flowerImplicitFunction() = default;


    // Member Functions

    virtual scalar value(const vector& p) const
    {
      scalar t=Foam::atan2(p.z(),p.x());
      scalar r=Foam::sqrt(Foam::sqr(p.z())+Foam::sqr(p.x()));
      return r+a*(Foam::cos(p*theta))//(-mag(project_ & (p - origin_)) + radius_)*scale_;
    }

    virtual vector grad(const vector& p) const
    {
        return -(project_ & (p - origin_))*scale_;
    }

    virtual scalar distanceToSurfaces(const vector& p) const
    {
        return mag(mag(project_ & (p - origin_)) - radius_)*scale_;
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace implicitFunctions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
