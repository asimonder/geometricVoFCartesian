/*---------------------------------------------------------------------------*\
    Part of interWYCenteredFoam (see the main file for the full
    description and license). Extracted for the conventional
    fragmented solver layout.
\*---------------------------------------------------------------------------*/

        // --- 2. viscous step (correction(dt) pattern) -----------------
        U += dtDim*gC;
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
            UEqn.solve();
            }
        }
        U -= dtDim*gC;
        U.correctBoundaryConditions();

