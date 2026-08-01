{
    // volume-fraction advection at the half-integer time levels
    // (staggered); no phase sources
    advector.advect();

    alpha2 = 1.0 - alpha1;
    rho == alpha1*rho1 + alpha2*rho2;

    if (turbulenceActive)
    {
        // mass-consistent face flux from the geometric alpha flux
        const surfaceScalarField& aPhi = advector.alphaPhi();
        rhoPhiPtr() = rho1*aPhi + rho2*(phi - aPhi);
    }

    // consistent momentum: the mixture momentum rides the same
    // geometric phase fluxes; sets Ustar and U = u*
    #include "mcvofAdvection.H"

    Info<< "Phase-1 volume fraction = "
        << alpha1.weightedAverage(mesh.Vsc()).value()
        << "  Min(" << alpha1.name() << ") = " << min(alpha1).value()
        << "  Max(" << alpha1.name() << ") = " << max(alpha1).value()
        << endl;
}
