/*---------------------------------------------------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     |                                                 |
|   \\  /    A nd           | www.openfoam.com                                |
|    \\/     M anipulation  |                                                 |
------------------------------------------------------------------------------
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

\*---------------------------------------------------------------------------*/

#include "interpolatedWaveImplicitFunction.H"
#include "addToRunTimeSelectionTable.H"
#include "Tuple2.H"
#include "fileName.H"
#include "tableBounds.H"
#include <vector>      // For std::vector
#include <fstream>     // For std::ifstream
#include <sstream>     // For std::stringstream
#include <algorithm>   // For std::sort

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace implicitFunctions
{
    defineTypeNameAndDebug(interpolatedWaveImplicitFunction, 0);
    addToRunTimeSelectionTable
    (
        implicitFunction,
        interpolatedWaveImplicitFunction,
        dict
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::implicitFunctions::interpolatedWaveImplicitFunction::interpolatedWaveImplicitFunction
(
    const dictionary& dict
)
:
    origin_(dict.get<vector>("origin")),
    direction_(normalised(dict.get<vector>("direction"))),
    up_(normalised(dict.get<vector>("up"))),
    h_(dict.get<scalar>("depth")),
    x0_(dict.get<scalar>("x0")),
    dataFile_(dict.get<fileName>("dataFile"))
{
    Info << "Constructing interpolatedWaveImplicitFunction..." << endl;
    buildEtaTable();
    Info << "Finished constructing interpolatedWaveImplicitFunction." << endl;
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void Foam::implicitFunctions::interpolatedWaveImplicitFunction::buildEtaTable()
{
    Info << "Building the x-eta interpolation table from file: " << dataFile_ << endl;

    std::ifstream file(dataFile_.c_str());
    if (!file.is_open())
    {
        FatalErrorInFunction
            << "Cannot open data file: " << dataFile_
            << abort(FatalError);
    }

    List<Tuple2<scalar, scalar>> tableData;
    std::string line;

    // Optional: Skip a header line if your CSV has one
    // std::getline(file, line);

    while (std::getline(file, line))
    {
        std::stringstream ss(line);
        std::string value;
        scalar x, eta;

        // Read x value from the first column
        if (std::getline(ss, value, ','))
        {
            // Use stod for C++11 and later for robust string-to-double conversion
            x = std::stod(value);
        }
        else
        {
            // Skip malformed lines
            continue;
        }

        // Read eta value from the second column
        if (std::getline(ss, value, ','))
        {
            eta = std::stod(value);
        }
        else
        {
            // Skip malformed lines
            continue;
        }
        
        // Add the point directly to our table data
        tableData.append(Tuple2<scalar, scalar>(x, eta));
    }
    file.close();

    Info << "  -> Read " << tableData.size() << " total points from the file." << endl;

    // --- Sort the data by x-value to ensure it's monotonic ---
    // This is a critical step for the interpolation table to work correctly.
    std::sort
    (
        tableData.begin(),
        tableData.end(),
        [](const Tuple2<scalar, scalar>& a, const Tuple2<scalar, scalar>& b)
        {
            return a.first() < b.first();
        }
    );

    // --- Create the OpenFOAM interpolation table ---
    etaTable_ = interpolationTable<scalar>
    (
        tableData,
        bounds::repeatableBounding::CLAMP, // Clamp values outside the table range
        fileName::null
    );

    Info << "  -> Successfully built interpolation table with "
         << etaTable_.size() << " points." << endl;
}

Foam::scalar Foam::implicitFunctions::interpolatedWaveImplicitFunction::value(const vector& p) const
{
    const scalar x = (p - origin_) & direction_;
    const scalar z = (p - origin_) & up_;

    scalar eta = etaTable_.interpolateValue((x-x0_)/h_);
    //    Info <<"x= "<<x-x0_<<", eta= "<<eta<< endl;
    
    return eta*h_ - z;
}

Foam::vector Foam::implicitFunctions::interpolatedWaveImplicitFunction::grad(const vector& p) const
{
    NotImplemented;
    return Zero;
}

Foam::scalar Foam::implicitFunctions::interpolatedWaveImplicitFunction::distanceToSurfaces(const vector& p) const
{
    NotImplemented;
    return 0;
}

// ************************************************************************* //


