/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2026 Asim Onder
     \\/     M anipulation  |
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

#include "geometricAdvector.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(geometricAdvector, 0);
    defineRunTimeSelectionTable(geometricAdvector, components);
}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::geometricAdvector> Foam::geometricAdvector::New
(
    volScalarField& alpha1,
    const surfaceScalarField& phi,
    const volVectorField& U
)
{
    const word advectorType
    (
        alpha1.mesh().solverDict(alpha1.name())
            .getOrDefault<word>("advector", "isoAdvector")
    );

    Info<< "Selecting geometric advector: " << advectorType << endl;

    auto* ctorPtr = componentsConstructorTable(advectorType);

    if (!ctorPtr)
    {
        FatalErrorInFunction
            << "Unknown advector " << advectorType << nl
            << "Valid advectors: "
            << componentsConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<geometricAdvector>(ctorPtr(alpha1, phi, U));
}


// ************************************************************************* //
