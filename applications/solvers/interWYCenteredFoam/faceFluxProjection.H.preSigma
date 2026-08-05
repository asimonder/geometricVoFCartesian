/*---------------------------------------------------------------------------*\
    Part of interWYCenteredFoam (see the main file for the full
    description and license). Extracted for the conventional
    fragmented solver layout.
\*---------------------------------------------------------------------------*/

        // --- 3. face flux + acceleration ------------------------------
        const surfaceScalarField rhof(fvc::interpolate(rho));
        // well-balanced face gravity acceleration [m/s2]
        const surfaceScalarField af(-ghf*fvc::snGrad(rho)/rhof);

        phi = (fvc::interpolate(U) & mesh.Sf())
            + dtDim*af*mesh.magSf();

        // --- 4. single exact face projection --------------------------
        const surfaceScalarField alphaf("alphaf", dtDim/rhof);
        constrainPressure(p_rgh, U, phi, alphaf);
        while (pimple.correctNonOrthogonal())
        {
            fvScalarMatrix pEqn
            (
                fvm::laplacian(alphaf, p_rgh) == fvc::div(phi)
            );
            pEqn.setReference
            (
                pRefCell, getRefCellValue(p_rgh, pRefCell)
            );
            pEqn.solve(p_rgh.select(pimple.finalNonOrthogonalIter()));
            if (pimple.finalNonOrthogonalIter())
            {
                phi -= pEqn.flux();
            }
        }
        #include "continuityErrs.H"

