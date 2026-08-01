    MRF.correctBoundaryVelocity(U);

    // Per-term temporal centering, read each outer iteration.
    // theta = 0 reproduces the stock Euler-implicit treatment.
    // theta = 0.5 is the trapezoid; the single interval flux
    // multiplies both convection halves, so uniform flow is
    // preserved for any thetaConv (weights sum to one).
    const scalar thetaConv =
        pimple.dict().getOrDefault<scalar>("thetaConv", 0.0);
    const scalar thetaVisc =
        pimple.dict().getOrDefault<scalar>("thetaVisc", 0.0);

    // BDF2 mode: pair backward ddt (set in fvSchemes) with the
    // extrapolated flux 1.5*rhoPhi - 0.5*rhoPhi.oldTime(), which
    // reproduces the three-level density history exactly.
    const bool bdf2Flux =
        pimple.dict().getOrDefault<bool>("bdf2Flux", false);

    surfaceScalarField rhoPhiEff
    (
        IOobject
        (
            "rhoPhi",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        rhoPhi
    );

    if (bdf2Flux)
    {
        rhoPhiEff = 1.5*rhoPhi - 0.5*rhoPhi.oldTime();
    }

    // Old-time velocity under the name "U", unregistered, so every
    // explicit half resolves the same fvSchemes entries as the
    // implicit half (scheme lookup goes by field name).
    volVectorField Uold
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        U.oldTime()
    );

    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U)
      + (1.0 - thetaConv)*fvm::div(rhoPhiEff, U)
      + thetaConv*fvc::div(rhoPhiEff, Uold)
      + MRF.DDt(rho, U)
      + (1.0 - thetaVisc)*turbulence->divDevRhoReff(rho, U)
      + thetaVisc*(turbulence->divDevRhoReff(rho, Uold) & Uold)
     ==
        fvOptions(rho, U)
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
