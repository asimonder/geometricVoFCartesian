    // U already holds the advected state u* = q*/rho(alpha_new)
    // (mcvofAdvection.H); the viscous step acts on u* + dt g and the
    // shift is removed after the solve (correction(dt) pattern)
    U += dtDim*gC;
    if (turbulenceActive)
    {
        // stock stress operator; the laminar branch below is the
        // original path and stays digit-identical when no model is
        // selected
        volVectorField Uprov("Uprov", U);
        fvVectorMatrix UEqn
        (
            fvm::Sp(rho/dtDim, U)
          + turbulence->divDevRhoReff(rho, U)
         ==
            rho/dtDim*Uprov
        );
        UEqn -= fvOptions(rho, U);
        fvOptions.constrain(UEqn);
        UEqn.solve();
    }
    else
    {
        const surfaceScalarField muf(mixture.muf());
        if (gMax(muf.primitiveField()) > VSMALL)
        {
            volVectorField Uprov("Uprov", U);
            fvVectorMatrix UEqn
            (
                fvm::Sp(rho/dtDim, U)
              - fvm::laplacian(muf, U)
             ==
                rho/dtDim*Uprov
              + fvc::div(muf*(mesh.Sf() & fvc::interpolate(dev2(Foam::T(fvc::grad(U))))))
            );
            UEqn -= fvOptions(rho, U);
            fvOptions.constrain(UEqn);
            UEqn.solve();
        }
    }
    U -= dtDim*gC;
    U.correctBoundaryConditions();

    fvOptions.correct(U);
