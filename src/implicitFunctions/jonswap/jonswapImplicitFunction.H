/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2025 Asim Onder
-------------------------------------------------------------------------------
Class
    Foam::implicitFunctions::jonswapImplicitFunction

Description
    Deep Water JONSWAP spectrum with exact periodic domain sizing.
    
    Inputs:
    - Hs, Tp, gamma
    - nWavesX, nWavesY (Domain size in # of peak wavelengths)
    
    The code calculates Lx = nWavesX * Lp to ensure the spectral peak 
    aligns exactly with a discrete wavenumber mode.
\*---------------------------------------------------------------------------*/

#ifndef implicitFunction_jonswapImplicitFunction_H
#define implicitFunction_jonswapImplicitFunction_H

#include "implicitFunction.H"
#include "List.H"
#include "vector.H"

namespace Foam
{
namespace implicitFunctions
{

class jonswapImplicitFunction
:
    public implicitFunction
{
private:

    // Private Data
        const scalar Hs_;
        const scalar Tp_;
        const scalar gamma_;
        const scalar beta_;
        const scalar theta0_;
        
        const scalar nWavesX_;
        const scalar nWavesY_;
        
        const label Nx_;
        const label Ny_;
        const label seed_;

        const vector origin_;
        const vector direction_;
        const vector up_;
        vector yDir_;

        scalar Lx_;
        scalar Ly_;

        // Spectral Modes (Written to file)
        List<scalar> amp_;
        List<scalar> kx_;
        List<scalar> ky_;
        List<scalar> phase_;
        List<scalar> omega_; // Added for velocity reconstruction

    // Private Member Functions
        void generateSpectrum();
        scalar jonswapSpectrum(scalar w, scalar wp) const;
        scalar spreadingFunction(scalar theta) const;
        void writeModes() const; // New function to dump data

public:

    TypeName("JONSWAP");

    explicit jonswapImplicitFunction(const dictionary& dict);
    virtual ~jonswapImplicitFunction() = default;

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
