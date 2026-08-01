/*---------------------------------------------------------------------------*\
    Part of interWYCenteredFoam (see the main file for the full
    description and license). Extracted for the conventional
    fragmented solver layout.
\*---------------------------------------------------------------------------*/

        // --- one-time alpha staggering at the dt plateau --------------
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

