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
    Foam::implicitFunctions::solitaryWaveClamondFructusImplicitFunction

Description
    Implicit function for a solitary wave based on the high-fidelity
    semi-analytical model of Clamond & Fructus (C. R. Mecanique 331, 2003).
    The wave profile is pre-calculated from the model's parameters and
    stored in an interpolation table for efficient evaluation.

SourceFiles
    solitaryWaveClamondFructus.C

\*---------------------------------------------------------------------------*/

#ifndef implicitFunction_solitaryWaveClamondFructusImplicitFunction_H
#define implicitFunction_solitaryWaveClamondFructusImplicitFunction_H

#include "implicitFunction.H"
#include "interpolationTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace implicitFunctions
{

/*---------------------------------------------------------------------------*\
           Class solitaryWaveClamondFructusImplicitFunction Declaration
\*---------------------------------------------------------------------------*/

class solitaryWaveClamondFructusImplicitFunction
:
    public implicitFunction
{
private:

    // Private Data

        //- Coordinate system vectors
        const vector up_;
        const vector direction_;
        const vector origin_;

        //- Physical parameters read from dictionary
        const scalar a_;         // H: Wave amplitude
        const scalar h_;         // d: Water depth
        const scalar x0_;

        //- Model parameters read from dictionary
        const scalar kappa_;
        const scalar nuParam_;
        const scalar A_;

        //- The pre-computed x-eta map for fast interpolation
        interpolationTable<scalar> etaTable_;


    // Private Member Functions

        //- A helper function to build the interpolation table
        void buildEtaTable();


public:

    //- Runtime type information
    TypeName("solitaryWaveClamondFructus");


    // Constructors

        //- Construct from dictionary
        explicit solitaryWaveClamondFructusImplicitFunction(const dictionary& dict);


    //- Destructor
    virtual ~solitaryWaveClamondFructusImplicitFunction() = default;


    // Member Functions

        //- Return the implicit function value at a point
        virtual scalar value(const vector& p) const;

        //- Return the gradient of the implicit function
        virtual vector grad(const vector& p) const;

        //- Return the distance to the nearest surface
        virtual scalar distanceToSurfaces(const vector& p) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace implicitFunctions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

