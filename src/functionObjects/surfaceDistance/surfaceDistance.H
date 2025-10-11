/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014-2016 OpenFOAM Foundation
    Copyright (C) 2016-2020 OpenCFD Ltd.
    Copyright (C) 2025 Asim Önder.
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
    Foam::functionObjects::surfaceDistance

Group
    grpFieldFunctionObjects

Description
    Computes the reconstructed distance function (surfaceDistance) field from a
    given volume-fraction field (typically alpha.water).

    The surfaceDistance is a signed distance field to the reconstructed interface
    generated using geometricVoF reconstruction schemes. It can be used
    for extracting wave-following surfaces, applying interface-aware
    models, or post-processing near-interface quantities.

    \f[
        \psi(\mathbf{x}) = \mathrm{signed\ distance\ to\ interface}
    \f]

    where
    \vartable
        \psi         | Reconstructed distance function [m]
        \alpha       | Volume fraction field [-]
    \endvartable

Operands:
    \table
      Operand        | Type           | Location
      input          | volScalarField | $FOAM_CASE/<time>/<alphaName>
      output field   | volScalarField | $FOAM_CASE/<time>/surfaceDistance
    \endtable

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    surfaceDistance1
    {
        // Mandatory entries
        type        surfaceDistance;
        libs        (libsurfaceDistance.so);

        // Optional entries
        alpha       alpha.water;   // name of volume fraction field
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                        | Type | Req'd | Dflt
      type         | Type name: surfaceDistance                     | word |  yes  | -
      libs         | Library name containing surfaceDistance        | word |  yes  | -
      alpha        | Name of volume-fraction field      | word |   no  | alpha.water
    \endtable

    The inherited entries are elaborated in:
     - \link functionObject.H \endlink
     - \link fvMeshFunctionObject.H \endlink

Minimal example by using the \c postProcess utility:
    \verbatim
        postProcess -func surfaceDistance
    \endverbatim

See also
    - Foam::functionObject
    - Foam::functionObjects::fvMeshFunctionObject
    - Foam::reconstructedDistanceFunction

SourceFiles
    surfaceDistance.C
\*---------------------------------------------------------------------------*/
#ifndef functionObjects_surfaceDistance_H
#define functionObjects_surfaceDistance_H

//#include "fvMeshFunctionObject.H"
#include "fieldExpression.H"
#include "reconstructedDistanceFunction.H"
#include "reconstructionSchemes.H"
#include "zoneDistribute.H"

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                           Class surfaceDistance Declaration
\*---------------------------------------------------------------------------*/

class surfaceDistance
:
    //    public fvMeshFunctionObject
    public fieldExpression
{
    // Private Data

        //- Name of the volume-fraction field
        word alphaName_;


    // Private Member Functions

        //- Calculate the reconstructed distance function and return true if successful
        virtual bool calc();

  //- Execute the functionObject
  // virtual bool execute();

  //- Write functionObject results
  //virtual bool write();
  
public:

    //- Runtime type information
    TypeName("surfaceDistance");


    // Constructors

        //- Construct from Time and dictionary
        surfaceDistance
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- No copy construct
        surfaceDistance(const surfaceDistance&) = delete;

        //- No copy assignment
        void operator=(const surfaceDistance&) = delete;

  
  inline label ij(const label i, const label j, const label Nx)
  {
    return j*Nx + i;
  }

  //- Destructor
    virtual ~surfaceDistance() = default;
};

} // End namespace functionObjects
} // End namespace Foam

#endif

// ************************************************************************* //
