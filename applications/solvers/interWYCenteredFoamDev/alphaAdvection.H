/*---------------------------------------------------------------------------*\
    Part of interWYCenteredFoam (see the main file for the full
    description and license). Extracted for the conventional
    fragmented solver layout.
\*---------------------------------------------------------------------------*/

        // --- 1. consistent advection (vof event) ----------------------
        q1 == alpha1*rho1*U;
        q2 == (scalar(1) - alpha1)*rho2*U;

        advector.advect(q1, q2);

        alpha2 = 1.0 - alpha1;
        mixture.correct();
        rho == alpha1*rho1 + alpha2*rho2;

        Info<< "Phase-1 volume fraction = "
            << alpha1.weightedAverage(mesh.Vsc()).value()
            << "  Min(" << alpha1.name() << ") = " << min(alpha1).value()
            << "  Max(" << alpha1.name() << ") = " << max(alpha1).value()
            << endl;

