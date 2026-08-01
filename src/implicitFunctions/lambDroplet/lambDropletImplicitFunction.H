/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2026 Asim Onder
-------------------------------------------------------------------------------
Class
    Foam::implicitFunctions::lambDropletImplicitFunction

Description
    Spherical droplet perturbed by a single spherical-harmonic mode, for the
    Lamb oscillating-droplet benchmark (Liu, Tolle, Bothe & Maric, JCP 493
    (2023) 112426, Sec. 4.2.4).

    Interface:

        r(theta,phi) = R0 + a * Y(theta,phi)

    with two shapes selected by the azimuthal number m:

        m = 0   Y = P_n(cos theta)              zonal, has a symmetry axis
        m = n   Y = sin^n(theta) cos(n phi)     sectoral, no symmetry axis

    Lamb's frequency depends on n only and not on m,

        omega_n^2 = n(n+1)(n-1)(n+2) sigma / ([(n+1) rho_d + n rho_a] R0^3),

    so both shapes carry the identical analytical reference. The sectoral
    variant is the genuinely three-dimensional one.

    Sign convention follows jonswapImplicitFunction: positive inside the
    liquid.

    Dictionary entries:

        origin      (0.02 0.02 0.02);   droplet centre
        R0          0.01;               unperturbed radius
        amplitude   0.00025;            perturbation amplitude, same units as R0
        n           2;                  mode number      (optional, default 2)
        m           0;                  0 or n           (optional, default 0)
        axis        (0 0 1);            polar axis       (optional)
        refDir      (1 0 0);            phi = 0 direction (optional)

\*---------------------------------------------------------------------------*/

#ifndef implicitFunction_lambDropletImplicitFunction_H
#define implicitFunction_lambDropletImplicitFunction_H

#include "implicitFunction.H"
#include "vector.H"

namespace Foam
{
namespace implicitFunctions
{

class lambDropletImplicitFunction
:
    public implicitFunction
{
private:

    // Private Data

        const vector origin_;
        const scalar R0_;
        const scalar amplitude_;
        const label  n_;
        const label  m_;

        vector axis_;
        vector refDir_;
        vector yDir_;

    // Private Member Functions

        //- Legendre polynomial P_n(x) by upward recurrence
        static scalar legendreP(const label n, const scalar x);

        //- The angular shape function Y(theta,phi)
        scalar shape(const vector& d, const scalar r) const;

public:

    TypeName("lambDroplet");

    explicit lambDropletImplicitFunction(const dictionary& dict);
    virtual ~lambDropletImplicitFunction() = default;

    virtual scalar value(const vector& p) const;

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

} // End namespace implicitFunctions
} // End namespace Foam

#endif
