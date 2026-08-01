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

\*---------------------------------------------------------------------------*/
#include "solitaryWaveClamondFructusImplicitFunction.H"
#include "addToRunTimeSelectionTable.H"
#include "Tuple2.H"
#include "fileName.H"
#include "tableBounds.H" 
#include <cmath>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace implicitFunctions
{
    // --- FIX: Use the correct class name ---
    defineTypeNameAndDebug(solitaryWaveClamondFructusImplicitFunction, 0);
    addToRunTimeSelectionTable
    (
        implicitFunction,
        solitaryWaveClamondFructusImplicitFunction,
        dict
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// --- FIX: Use the correct class name ---
Foam::implicitFunctions::solitaryWaveClamondFructusImplicitFunction::solitaryWaveClamondFructusImplicitFunction
(
    const dictionary& dict
)
:
    // ... (initializer list is the same, but ensure order matches your .H file) ...
    origin_(dict.get<vector>("origin")),
    direction_(normalised(dict.get<vector>("direction"))),
    up_(normalised(dict.get<vector>("up"))),
    x0_(dict.get<scalar>("x0")),
    a_(dict.get<scalar>("amplitude")),
    h_(dict.get<scalar>("depth")),
    kappa_(dict.get<scalar>("kappa")),
    nuParam_(dict.get<scalar>("nuParam")),
    A_(dict.get<scalar>("A"))
{
    Info << "Constructing solitaryWaveClamondFructusImplicitFunction..." << endl;
    buildEtaTable();
    Info << "Finished constructing solitaryWaveClamondFructusImplicitFunction." << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void Foam::implicitFunctions::solitaryWaveClamondFructusImplicitFunction::buildEtaTable()
{
    Info << "Building the x-eta interpolation table..." << endl;

    const int nPoints = 801;
    List<scalar> theta_half(nPoints);
    List<scalar> eta_half(nPoints);

    const scalar chi = kappa_ * (1.0/A_) * Foam::sin(nuParam_);

    for (int i = 0; i < nPoints; ++i)
    {
        scalar eta = a_ * Foam::pow(static_cast<scalar>(i) / (nPoints - 1), 2.0);
        eta_half[i] = eta;

        scalar calculated_theta = -1.0; // Default to invalid

        // Handle the crest point as a special case
        if (mag(eta - a_) < 1e-9)
        {
            calculated_theta = 0.0;
        }
        else
        {
            const scalar arg_coth = chi * eta;
            if (mag(arg_coth) > 1e-9)
            {
                const scalar rhs_term1 = Foam::sin(nuParam_) * Foam::sin(kappa_ * (eta + h_)) / Foam::tanh(arg_coth);
                const scalar rhs_term2 = Foam::cos(nuParam_) * Foam::cos(kappa_ * (eta + h_));
                const scalar rhs = rhs_term1 - rhs_term2;

                if (rhs >= 1.0)
                {
                    calculated_theta = std::acosh(rhs) / kappa_;
                }

            }
        }
        theta_half[i] = calculated_theta;
    }

    List<Tuple2<scalar, scalar>> tableData;

    for (int i = nPoints - 1; i > 0; --i)
    {
        if (theta_half[i] >= 0)
        {
            tableData.append(Tuple2<scalar, scalar>(-theta_half[i], eta_half[i]));
        }
    }

    for (int i = 0; i < nPoints; ++i)
    {
        if (theta_half[i] >= 0)
        {
            tableData.append(Tuple2<scalar, scalar>(theta_half[i], eta_half[i]));
        }
    }

    std::sort
    (
        tableData.begin(),
        tableData.end(),
        [](const Tuple2<scalar, scalar>& a, const Tuple2<scalar, scalar>& b)
        {
            return a.first() < b.first();
        }
    );
	
    // Use the correct constructor with the bounds::repeatableBounding enum.
    etaTable_ = interpolationTable<scalar>
    (
        tableData,
        bounds::repeatableBounding::CLAMP,
        fileName::null
    );

    Info << "  -> Successfully built interpolation table with "
         << etaTable_.size() << " points." << endl;
    Info<<etaTable_<<endl;
}


Foam::scalar Foam::implicitFunctions::solitaryWaveClamondFructusImplicitFunction::value(const vector& p) const
{
    const scalar x = (p - origin_) & direction_;
    const scalar z = (p - origin_) & up_;

    scalar eta = etaTable_.interpolateValue(x-x0_);
    Info <<"x= "<<x-x0_<<", eta= "<<eta<< endl;
    
    return eta - z;
}


// --- FIX: Use the correct class name ---
Foam::vector Foam::implicitFunctions::solitaryWaveClamondFructusImplicitFunction::grad(const vector& p) const
{
    NotImplemented;
    return Zero;
}


// --- FIX: Use the correct class name ---
Foam::scalar Foam::implicitFunctions::solitaryWaveClamondFructusImplicitFunction::distanceToSurfaces(const vector& p) const
{
    NotImplemented;
    return 0;
}


// ************************************************************************* //
