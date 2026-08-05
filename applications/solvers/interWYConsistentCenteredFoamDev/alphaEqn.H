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

    // consistent advection: phase momenta ride the alpha fluxes
    Upre = U;

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

    Info<< "Phase-1 volume fraction = "
        << alpha1.weightedAverage(mesh.Vsc()).value()
        << "  Min(" << alpha1.name() << ") = " << min(alpha1).value()
        << "  Max(" << alpha1.name() << ") = " << max(alpha1).value()
        << endl;
}
