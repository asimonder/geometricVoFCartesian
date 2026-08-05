    // consistent momentum recovery: momentum is advected on every
    // cell by the sweeps (geometric on the lattice, per-phase
    // Lax-Wendroff in the pure-phase bulk), so the one-fluid
    // velocity is recovered everywhere from the phase momenta
    // DevPhiAdv: pre-viscous state, so the face predictor can
    // carry the viscous increment without the full cell field
    volVectorField UpreVisc
    (
        IOobject("UpreVisc", runTime.timeName(), mesh,
                 IOobject::NO_READ, IOobject::NO_WRITE,
                 IOobject::NO_REGISTER),
        U
    );

    U.primitiveFieldRef() =
        (q1.primitiveField() + q2.primitiveField())
       /rho.primitiveField();
    U.correctBoundaryConditions();

    // implicit viscous step on u + dt g (correction(dt) pattern)
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
            // fvOptions(rho, U) subtracted into the matrix; the reference's
            // "== fvOptions(rho, U)" expands to LHS - fvOptions via operator-=
            UEqn -= fvOptions(rho, U);
            fvOptions.constrain(UEqn);
            UEqn.solve();
        }
    }
    U -= dtDim*gC;
    U.correctBoundaryConditions();

    // predictor-stage fvOptions correction (reference UEqn.H);
    // the post-projection correction is in pEqn.H
    fvOptions.correct(U);
