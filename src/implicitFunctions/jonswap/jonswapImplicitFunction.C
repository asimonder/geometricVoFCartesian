/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2025 Asim Onder
-------------------------------------------------------------------------------
\*---------------------------------------------------------------------------*/

#include "jonswapImplicitFunction.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"
#include "Random.H"
#include "OFstream.H"
#include "OSspecific.H"

using namespace Foam::constant;

namespace Foam
{
namespace implicitFunctions
{

defineTypeNameAndDebug(jonswapImplicitFunction, 0);
addToRunTimeSelectionTable(implicitFunction, jonswapImplicitFunction, dict);

jonswapImplicitFunction::jonswapImplicitFunction(const dictionary& dict)
:
    Hs_(dict.get<scalar>("Hs")),
    Tp_(dict.get<scalar>("Tp")),
    gamma_(dict.getOrDefault<scalar>("gamma", 3.3)),
    beta_(dict.getOrDefault<scalar>("beta", 0.174532925)),
    theta0_(dict.getOrDefault<scalar>("theta0", 0.0)),
    nWavesX_(dict.get<scalar>("nWavesX")),
    nWavesY_(dict.get<scalar>("nWavesY")),
    Nx_(dict.get<label>("Nx")),
    Ny_(dict.get<label>("Ny")),
    seed_(dict.getOrDefault<label>("seed", 12345)),
    origin_(dict.get<vector>("origin")),
    direction_(dict.get<vector>("direction").normalise()),
    up_(dict.get<vector>("up").normalise()),
    yDir_(up_ ^ direction_)
{
    generateSpectrum();
    writeModes(); // Dump to file immediately
}

void jonswapImplicitFunction::generateSpectrum()
{
    const scalar pi = mathematical::pi;
    const scalar g = 9.81; 

    const scalar wp = 2.0*pi/Tp_;
    const scalar kp = (wp*wp)/g; 
    const scalar Lp = 2.0*pi/kp;

    Lx_ = nWavesX_ * Lp;
    Ly_ = nWavesY_ * Lp;

    Info<< "Initializing JONSWAP (Hs=" << Hs_ << ", Tp=" << Tp_ << ")" << endl;

    scalar targetVariance = pow(Hs_/4.0, 2);
    scalar scaleFactor = targetVariance; 
    
    scalar error = 1.0;
    label iter = 0;
    const label maxIter = 20;
    const scalar tol = 1e-4;

    // Phase 1: Scaling Loop
    while (error > tol && iter < maxIter)
    {
        scalar totalVariance = 0.0;
        for (label i = 1; i <= Nx_; ++i)
        {
            scalar kx_val = i * (2.0*pi/Lx_);
            for (label j = -Ny_/2; j < Ny_/2; ++j)
            {
                scalar ky_val = j * (2.0*pi/Ly_);
                scalar k = sqrt(kx_val*kx_val + ky_val*ky_val);
                scalar theta = atan2(ky_val, kx_val);
                scalar w = sqrt(g * k);
                scalar Cg = 0.5 * (w/k);

                scalar Sw = jonswapSpectrum(w, wp) * scaleFactor; 
                scalar Dtheta = spreadingFunction(theta);
                scalar dkx = 2.0*pi/Lx_;
                scalar dky = 2.0*pi/Ly_;
                
                scalar E_bin = (Sw * Dtheta) * (Cg/k) * dkx * dky;
                if (E_bin > 1e-20) totalVariance += E_bin;
            }
        }

        if (totalVariance > 1e-15) {
            scalar ratio = targetVariance / totalVariance;
            scaleFactor *= ratio;
            error = mag(1.0 - ratio);
        } else { scaleFactor *= 2.0; }
        iter++;
    }

    // Phase 2: Generation
    Random rndGen(seed_);
    DynamicList<scalar> tmpAmp, tmpKx, tmpKy, tmpPhase, tmpOmega;
    
    for (label i = 1; i <= Nx_; ++i)
    {
        scalar kx_val = i * (2.0*pi/Lx_);
        for (label j = -Ny_/2; j < Ny_/2; ++j)
        {
            scalar ky_val = j * (2.0*pi/Ly_);
            scalar k = sqrt(kx_val*kx_val + ky_val*ky_val);
            scalar theta = atan2(ky_val, kx_val);
            scalar w = sqrt(g * k);
            scalar Cg = 0.5 * (w/k);

            scalar Sw = jonswapSpectrum(w, wp) * scaleFactor; 
            scalar Dtheta = spreadingFunction(theta);
            scalar dkx = 2.0*pi/Lx_;
            scalar dky = 2.0*pi/Ly_;
            scalar E_bin = (Sw * Dtheta) * (Cg/k) * dkx * dky;

            if (E_bin > 1e-20)
            {
                tmpAmp.append(sqrt(2.0 * E_bin));
                tmpKx.append(kx_val);
                tmpKy.append(ky_val);
                tmpOmega.append(w);
                tmpPhase.append(rndGen.sample01<scalar>() * 2.0 * pi);
            }
        }
    }

    amp_.transfer(tmpAmp);
    kx_.transfer(tmpKx);
    ky_.transfer(tmpKy);
    phase_.transfer(tmpPhase);
    omega_.transfer(tmpOmega);
    
    Info<< "  Generated " << amp_.size() << " modes." << endl;
}

void jonswapImplicitFunction::writeModes() const
{
    fileName path("constant/waveModes.dat");
    OFstream os(path);
    
    if (!os.good())
    {
        FatalErrorInFunction << "Cannot open " << path << " for writing." << exit(FatalError);
    }

    Info<< "Writing wave modes to " << path << "..." << endl;

    // Write total number of scalars (5 * nModes)
    os << amp_.size() * 5 << nl;
    os << "(" << nl;
    
    forAll(amp_, i)
    {
        // Write raw numbers, space separated, no parentheses
        os << amp_[i] << " " 
           << kx_[i] << " " 
           << ky_[i] << " " 
           << omega_[i] << " " 
           << phase_[i] << nl;
    }
    os << ");" << nl;
}

scalar jonswapImplicitFunction::jonswapSpectrum(scalar w, scalar wp) const
{
    scalar sigma = (w <= wp) ? 0.07 : 0.09;
    scalar arg = -0.5 * pow((w - wp) / (sigma * wp), 2);
    scalar Spm = pow(w, -5) * exp(-1.25 * pow(wp/w, 4));
    return Spm * pow(gamma_, exp(arg));
}

scalar jonswapImplicitFunction::spreadingFunction(scalar theta) const
{
    const scalar pi = mathematical::pi;
    scalar diff = theta - theta0_;
    while (diff <= -pi) diff += 2*pi;
    while (diff > pi) diff -= 2*pi;
    
    if (mag(diff) <= beta_)
    {
        return (1.0/beta_) * pow(cos(pi*diff/(2.0*beta_)), 2);
    }
    return 0.0;
}

scalar jonswapImplicitFunction::value(const vector& p) const
{
    vector d = p - origin_;
    scalar x = d & direction_;
    scalar y = d & yDir_;
    scalar z = d & up_;

    scalar eta = 0.0;
    const label nModes = amp_.size();
    const scalar* ampPtr = amp_.begin();
    const scalar* kxPtr = kx_.begin();
    const scalar* kyPtr = ky_.begin();
    const scalar* phPtr = phase_.begin();

    for (label i = 0; i < nModes; ++i)
    {
        eta += ampPtr[i] * cos(kxPtr[i]*x + kyPtr[i]*y + phPtr[i]);
    }

    return eta - z;
}

} // End namespace implicitFunctions
} // End namespace Foam

