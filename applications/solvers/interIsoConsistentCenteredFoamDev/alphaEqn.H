{
    // volume-fraction advection at the half-integer time levels
    // (staggered); no phase sources
    #include "alphaSuSp.H"

    // staggered-interval maintenance. Alpha exits each step half a
    // step ahead of the new time level. That offset is exact only at
    // constant step; when the step changes, keeping it requires the
    // advection interval (3 dt_new - dt_old)/2, which equals dt at
    // constant step. Implemented by scaling phi during the advection.
    // The scale is clamped to [0.5, 1.5] so the advector's Courant
    // bound survives an abrupt step collapse; the residual offset of
    // a clamped event is absorbed by later steps.
    scalar intervalScale = 1.0;
    if (staggerAlpha && staggerDone)
    {
        const scalar dtNew = runTime.deltaTValue();
        const scalar dtOld = runTime.deltaT0Value();
        if (dtOld > VSMALL && dtNew > VSMALL)
        {
            intervalScale = min(max(0.5*(3.0 - dtOld/dtNew), 0.5), 1.5);
        }
    }
    const bool intervalScaled = mag(intervalScale - 1.0) > 1e-12;
    if (intervalScaled) { phi *= intervalScale; }

    Upre = U;

    advector.advect(Sp, Su);

    alpha2 = 1.0 - alpha1;
    rho == alpha1*rho1 + alpha2*rho2;

    if (turbulenceActive)
    {
        // mass-consistent face flux from the geometric alpha flux
        const surfaceScalarField& aPhi = advector.alphaPhi();
        rhoPhiPtr() = rho1*aPhi + rho2*(phi - aPhi);
    }

    // consistent momentum: the mixture momentum rides the same
    // geometric phase fluxes; sets Ustar and U = u*. Runs on the
    // scaled phi so the gas flux phi - alphaPhi stays coherent with
    // the scaled liquid flux.
    #include "mcvofAdvection.H"

    if (intervalScaled) { phi /= intervalScale; }

    Info<< "Phase-1 volume fraction = "
        << alpha1.weightedAverage(mesh.Vsc()).value()
        << "  Min(" << alpha1.name() << ") = " << min(alpha1).value()
        << "  Max(" << alpha1.name() << ") = " << max(alpha1).value()
        << endl;
}
