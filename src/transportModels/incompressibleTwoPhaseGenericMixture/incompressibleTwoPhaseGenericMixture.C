/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "incompressibleTwoPhaseGenericMixture.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(incompressibleTwoPhaseGenericMixture, 0);
}


// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

void Foam::incompressibleTwoPhaseGenericMixture::calcNu()
{
    nuModel1_->correct();
    nuModel2_->correct();

    const volScalarField limitedAlpha1
    (
        "limitedAlpha1",
        clamp(alpha1_, zero_one{})
    );

    // Average kinematic viscosity calculated from dynamic viscosity
    nu_ = mu()/(limitedAlpha1*rho1_ + (scalar(1) - limitedAlpha1)*rho2_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::incompressibleTwoPhaseGenericMixture::incompressibleTwoPhaseGenericMixture
(
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    IOdictionary
    (
        IOobject
        (
            "transportProperties",
            U.time().constant(),
            U.db(),
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    twoPhaseMixture(U.mesh(), *this),

    nuModel1_
    (
        viscosityModel::New
        (
            "nu1",
            subDict(phase1Name_),
            U,
            phi
        )
    ),
    nuModel2_
    (
        viscosityModel::New
        (
            "nu2",
            subDict(phase2Name_),
            U,
            phi
        )
    ),

    rho1_("rho", dimDensity, nuModel1_->viscosityProperties()),
    rho2_("rho", dimDensity, nuModel2_->viscosityProperties()),

    U_(U),
    phi_(phi),

    nu_
    (
        IOobject
        (
            "nu",
            U_.time().timeName(),
            U_.db()
        ),
        U_.mesh(),
        dimensionedScalar(dimViscosity, Zero)
    )
{
  phaseInterpolation_ = this->lookupOrDefault<word>("phaseInterpolation", "linear");
  if (phaseInterpolation_ != "linear" && phaseInterpolation_ != "harmonic")
    {
      FatalIOErrorInFunction(*this)
        << "Unknown phaseInterpolation type: " << phaseInterpolation_ << nl
        << "Valid options are:" << nl
        << "    linear"  << nl
        << "    harmonic" << nl
        << exit(FatalIOError);
    }
  
  calcNu();
  
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::incompressibleTwoPhaseGenericMixture::mu() const
{
    const volScalarField limitedAlpha1
    (
        clamp(alpha1_, zero_one{})
    );

    const volScalarField mu1 = (rho1_*nuModel1_->nu())();
    const volScalarField mu2 = (rho2_*nuModel2_->nu())();

    if (phaseInterpolation_ == "harmonic")
    {
        return volScalarField::New
        (
            "mu",
            IOobject::NO_REGISTER,
            (mu1*mu2)/(limitedAlpha1*mu2 + (scalar(1) - limitedAlpha1)*mu1)
        );
    }
    else
    {
        return volScalarField::New
        (
            "mu",
            IOobject::NO_REGISTER,
            limitedAlpha1*mu1 + (scalar(1) - limitedAlpha1)*mu2
        );
    }
}


Foam::tmp<Foam::scalarField>
Foam::incompressibleTwoPhaseGenericMixture::mu(const label patchI) const
{

    return mu()().boundaryField()[patchI];
}

Foam::tmp<Foam::surfaceScalarField>
Foam::incompressibleTwoPhaseGenericMixture::muf() const
{
    const surfaceScalarField alpha1f
    (
        clamp(fvc::interpolate(alpha1_), zero_one{})
    );

    const surfaceScalarField mu1f = (rho1_ * fvc::interpolate(nuModel1_->nu()))();
    const surfaceScalarField mu2f = (rho2_ * fvc::interpolate(nuModel2_->nu()))();

    if (phaseInterpolation_ == "harmonic")
    {
        return surfaceScalarField::New
        (
            "muf",
            IOobject::NO_REGISTER,
            (mu1f*mu2f) / (alpha1f*mu2f + (scalar(1) - alpha1f)*mu1f)
        );
    }
    else
    {
        return surfaceScalarField::New
        (
            "muf",
            IOobject::NO_REGISTER,
            alpha1f*mu1f + (scalar(1) - alpha1f)*mu2f
        );
    }
}

Foam::tmp<Foam::surfaceScalarField>
Foam::incompressibleTwoPhaseGenericMixture::nuf() const
{
    const surfaceScalarField alpha1f
    (
        clamp(fvc::interpolate(alpha1_), zero_one{})
    );

    const surfaceScalarField mu1f = (rho1_ * fvc::interpolate(nuModel1_->nu()))();
    const surfaceScalarField mu2f = (rho2_ * fvc::interpolate(nuModel2_->nu()))();
    const surfaceScalarField rhof =
      (alpha1f*rho1_ + (scalar(1) - alpha1f)*rho2_)();

    if (phaseInterpolation_ == "harmonic")
    {
        return surfaceScalarField::New
        (
            "nuf",
            IOobject::NO_REGISTER,
            (mu1f*mu2f) / (alpha1f*mu2f + (scalar(1) - alpha1f)*mu1f) / rhof
        );
    }
    else
    {
        return surfaceScalarField::New
        (
            "nuf",
            IOobject::NO_REGISTER,
            (alpha1f*mu1f + (scalar(1) - alpha1f)*mu2f) / rhof
        );
    }
}

bool Foam::incompressibleTwoPhaseGenericMixture::read()
{
    if (regIOobject::read())
    {
        if
        (
            nuModel1_().read
            (
                subDict(phase1Name_ == "1" ? "phase1": phase1Name_)
            )
         && nuModel2_().read
            (
                subDict(phase2Name_ == "2" ? "phase2": phase2Name_)
            )
        )
        {
            nuModel1_->viscosityProperties().readEntry("rho", rho1_);
            nuModel2_->viscosityProperties().readEntry("rho", rho2_);

            return true;
        }
    }

    return false;
}


// ************************************************************************* //
