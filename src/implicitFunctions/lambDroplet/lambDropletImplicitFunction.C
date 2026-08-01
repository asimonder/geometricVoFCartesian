/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2026 Asim Onder
-------------------------------------------------------------------------------
\*---------------------------------------------------------------------------*/

#include "lambDropletImplicitFunction.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"

namespace Foam
{
namespace implicitFunctions
{

defineTypeNameAndDebug(lambDropletImplicitFunction, 0);
addToRunTimeSelectionTable(implicitFunction, lambDropletImplicitFunction, dict);

lambDropletImplicitFunction::lambDropletImplicitFunction(const dictionary& dict)
:
    origin_(dict.get<vector>("origin")),
    R0_(dict.get<scalar>("R0")),
    amplitude_(dict.get<scalar>("amplitude")),
    n_(dict.getOrDefault<label>("n", 2)),
    m_(dict.getOrDefault<label>("m", 0)),
    axis_(dict.getOrDefault<vector>("axis", vector(0, 0, 1))),
    refDir_(dict.getOrDefault<vector>("refDir", vector(1, 0, 0))),
    yDir_(Zero)
{
    if (n_ < 2)
    {
        FatalErrorInFunction
            << "Mode number n must be at least 2, got " << n_
            << exit(FatalError);
    }

    if (m_ != 0 && m_ != n_)
    {
        FatalErrorInFunction
            << "Only the zonal mode m = 0 and the sectoral mode m = n are"
            << " implemented. Got n = " << n_ << ", m = " << m_
            << exit(FatalError);
    }

    axis_.normalise();

    // Gram-Schmidt: make refDir orthogonal to the axis, then complete the triad
    refDir_ -= (refDir_ & axis_)*axis_;

    if (mag(refDir_) < SMALL)
    {
        FatalErrorInFunction
            << "refDir " << dict.getOrDefault<vector>("refDir", vector(1, 0, 0))
            << " is parallel to axis " << axis_
            << exit(FatalError);
    }

    refDir_.normalise();
    yDir_ = axis_ ^ refDir_;

    Info<< "lambDroplet: R0 = " << R0_
        << ", amplitude = " << amplitude_
        << ", n = " << n_ << ", m = " << m_
        << (m_ == 0 ? " (zonal, axisymmetric)" : " (sectoral, fully 3D)")
        << ", origin = " << origin_ << endl;
}


scalar lambDropletImplicitFunction::legendreP(const label n, const scalar x)
{
    if (n == 0)
    {
        return 1;
    }

    scalar pkm1 = 1;    // P_0
    scalar pk = x;      // P_1

    for (label k = 1; k < n; ++k)
    {
        const scalar pkp1 = ((2*k + 1)*x*pk - k*pkm1)/(k + 1);
        pkm1 = pk;
        pk = pkp1;
    }

    return pk;
}


scalar lambDropletImplicitFunction::shape(const vector& d, const scalar r) const
{
    const scalar cosTheta = min(1.0, max(-1.0, (d & axis_)/r));

    if (m_ == 0)
    {
        return legendreP(n_, cosTheta);
    }

    // sectoral: sin^n(theta) cos(n phi)
    const scalar sinTheta = sqrt(max(0.0, 1 - cosTheta*cosTheta));
    const vector dPerp = d - (d & axis_)*axis_;
    const scalar phi = atan2(dPerp & yDir_, dPerp & refDir_);

    return pow(sinTheta, scalar(n_))*cos(n_*phi);
}


scalar lambDropletImplicitFunction::value(const vector& p) const
{
    const vector d = p - origin_;
    const scalar r = mag(d);

    if (r < SMALL)
    {
        // at the centre the angular part is undefined; deep inside either way
        return R0_;
    }

    return (R0_ + amplitude_*shape(d, r)) - r;
}


} // End namespace implicitFunctions
} // End namespace Foam
