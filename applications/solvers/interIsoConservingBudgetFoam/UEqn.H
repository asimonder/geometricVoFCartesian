    MRF.correctBoundaryVelocity(U);

    // MCVOF: momentum advection was applied explicitly (U holds u*).
    // Remaining implicit balance from the advected state:
    //   rho^{n+1} (U - u*)/dt = viscous + forces - grad(p)
    // The convection term does NOT appear here by construction.
    fvVectorMatrix UEqn
    (
        fvm::Sp(rho/runTime.deltaT(), U)
      + MRF.DDt(rho, U)
      + turbulence->divDevRhoReff(rho, U)
     ==
        rho*Ustar/runTime.deltaT()
      + fvOptions(rho, U)
    );

    UEqn.relax();

    fvOptions.constrain(UEqn);

    if (pimple.momentumPredictor())
    {
        solve
        (
            UEqn
	    ==
            fvc::reconstruct
            (
                (
		 interfaceForce.surfaceTensionForce()
		 - ghf*fvc::snGrad(rho)
                 - fvc::snGrad(p_rgh)
                ) * mesh.magSf()
            )
        );

        fvOptions.correct(U);
    }
