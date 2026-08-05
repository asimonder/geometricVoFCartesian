{
    // one-time dt/2 alpha staggering at the dt plateau (staggered-in-time
    // volume fraction of the centered cycle)
    if (staggerAlpha && !staggerDone)
    {
        const scalar dtNow = runTime.deltaTValue();
        if (dtPrev > 0 && mag(dtNow - dtPrev) < 1e-3*dtNow)
        {
            Info<< "wy: dt plateau; staggering alpha by dt/2 = "
                << 0.5*dtNow << " s" << nl;
            phi *= 0.5;
            advector.advect();
            phi *= 2.0;
            alpha2 = 1.0 - alpha1;
            alpha2.correctBoundaryConditions();
            mixture.correct();
            rho == alpha1*rho1 + alpha2*rho2;
            staggerDone = true;
        }
        dtPrev = dtNow;
    }

    // staggered-interval maintenance. Alpha exits each step half a
    // step ahead of the new time level. That offset is exact only at
    // constant step; when the step changes, keeping it requires the
    // advection interval (3 dt_new - dt_old)/2, which equals dt at
    // constant step. Implemented by scaling phi during the advection,
    // which scales the sweep volumes, the compression term, and the
    // momentum transport coherently. The scale is clamped to
    // [0.5, 1.5] so the per-sweep Courant bound survives an abrupt
    // step collapse; the residual offset of a clamped event is
    // absorbed by later steps.
    scalar intervalScale = 1.0;
    if (false)  // NoSI experiment: interval maintenance disabled
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

    // consistent advection: phase momenta ride the alpha fluxes
    q1 == alpha1*rho1*U;
    q2 == (1.0 - alpha1)*rho2*U;

    advector.advect(q1, q2);

    alpha2 = 1.0 - alpha1;
    rho == alpha1*rho1 + alpha2*rho2;

    if (turbulenceActive)
    {
        // mass-consistent face flux from the geometric alpha flux
        const surfaceScalarField& aPhi = advector.alphaPhi();
        rhoPhiPtr() = rho1*aPhi + rho2*(phi - aPhi);
    }

    if (intervalScaled) { phi /= intervalScale; }

    Info<< "Phase-1 volume fraction = "
        << alpha1.weightedAverage(mesh.Vsc()).value()
        << "  Min(" << alpha1.name() << ") = " << min(alpha1).value()
        << "  Max(" << alpha1.name() << ") = " << max(alpha1).value()
        << endl;
}
