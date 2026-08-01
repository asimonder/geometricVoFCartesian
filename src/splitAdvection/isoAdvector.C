/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2026 Asim Onder
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM (GPL v3 or later; see
    <http://www.gnu.org/licenses/>).

\*---------------------------------------------------------------------------*/

#include "isoAdvector.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
    defineTypeNameAndDebug(isoAdvector, 0);
    addToRunTimeSelectionTable
    (
        geometricAdvector,
        isoAdvector,
        components
    );
}

// ************************************************************************* //
