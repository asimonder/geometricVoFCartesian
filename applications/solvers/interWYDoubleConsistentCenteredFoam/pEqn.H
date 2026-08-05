{
    // Double approximate projection (Basilisk
    // navier-stokes/double-projection.h, Almgren et al. 2000). Two
    // pressures: pAux projects the transport flux and is the only
    // pressure the face system feels, so the transport projection no
    // longer sees the divergence history of the cell field; p_rgh is
    // obtained by projecting the per-step UPDATE system Af and
    // supplies the centered gradient gC applied at the NEXT step.

    // HARMONIC face density: the weight family of the droplet-passing
    // collocated solvers (rAUf = interp(dt/rho) in interIsoCartFoam)
    const surfaceScalarField rhof(1.0/fvc::interpolate(1.0/rho));
    const surfaceScalarField rAUf("rAUf", dtDim/rhof);

    // balanced face acceleration: gravity + surface tension share rhof
    // and the discrete operators with the pressure gradient
    surfaceScalarField af(-ghf*fvc::snGrad(rho)/rhof);
    if (stActive)
    {
        af += interfaceForcePtr->surfaceTensionForce()/rhof;
    }
    if (airForcing)
    {
        // air-partial force density at faces, same weights as gravity
        // and surface tension; the cell correction gC returns it to
        // cells. Force density (1-alpha_f)*G is the air-side stress
        // budget on every face including mixed ones (the Basilisk
        // original adds a bare kinematic acceleration instead and
        // over-delivers at mixed faces by the local density ratio).
        af += (1.0 - fvc::interpolate(alpha1))
             *((airForce & mesh.Sf())/mesh.magSf())/rhof;
    }

    // the per-step update system: advection + viscous face increment
    // plus the face forces (double-projection step 1). U here holds
    // the advected and diffused state, U.oldTime() the corrected
    // velocity at the start of the step.
    const surfaceScalarField Af
    (
        "Af",
        (fvc::interpolate(U - U.oldTime()) & mesh.Sf())
      + dtDim*af*mesh.magSf()
    );

    // permanent centered correction with the PREVIOUS step's gradient
    // (the correction(dt) that double-projection.h makes permanent)
    U += dtDim*gC;
    U.correctBoundaryConditions();

    // ---- transport projection with the auxiliary pressure ---------
    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        fvc::interpolate(U) & mesh.Sf()
    );

    constrainPressure(pAux, U, phiHbyA, rAUf);

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pAuxEqn
        (
            fvm::laplacian(rAUf, pAux) == fvc::div(phiHbyA)
        );
        pAuxEqn.setReference
        (
            pRefCell, getRefCellValue(pAux, pRefCell)
        );
        pAuxEqn.solve(mesh.solverDict(p_rgh.name()));
        if (pimple.finalNonOrthogonalIter())
        {
            phi = phiHbyA - pAuxEqn.flux();
        }
    }
    #include "continuityErrs.H"

    // approximate cell correction with grad(pAux) alone
    // (double-projection step 4): same gather pairing as production,
    // scratch target so gC keeps its own role
    {
        const surfaceScalarField hf1(-rhof*(rAUf/dtDim)*fvc::snGrad(pAux));
        faceToCellAvg(hf1, gScratch);
        {
            vectorField& gcf = gScratch.primitiveFieldRef();
            const scalarField& rhoc = rho.primitiveField();
            forAll(gcf, celli)
            {
                gcf[celli] /= rhoc[celli];
            }
        }
        if (zonal)
        {
            const volVectorField gRec(fvc::reconstruct(hf1*mesh.magSf()));
            const boolList& zc = splitPtr->inZone();
            vectorField& gcf = gScratch.primitiveFieldRef();
            const scalarField& rhoc = rho.primitiveField();
            forAll(gcf, celli)
            {
                if (!zc[celli]) gcf[celli] = gRec[celli]/rhoc[celli];
            }
        }
        gScratch.correctBoundaryConditions();

        U += dtDim*gScratch;
        U.correctBoundaryConditions();
    }

    fvOptions.correct(U);

    // ---- second projection: the update system defines the true
    // pressure and the centered gradient for the next step ----------
    constrainPressure(p_rgh, U, Af, rAUf);

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix p_rghEqn
        (
            fvm::laplacian(rAUf, p_rgh) == fvc::div(Af)
        );
        p_rghEqn.setReference
        (
            pRefCell, getRefCellValue(p_rgh, pRefCell)
        );
        p_rghEqn.solve(p_rgh.select(pimple.finalNonOrthogonalIter()));
    }

    const surfaceScalarField gf(af - (rAUf/dtDim)*fvc::snGrad(p_rgh));
    const surfaceScalarField hf(rhof*gf);
    faceToCellAvg(hf, gC);
    {
        vectorField& gcf = gC.primitiveFieldRef();
        const scalarField& rhoc = rho.primitiveField();
        forAll(gcf, celli)
        {
            gcf[celli] /= rhoc[celli];
        }
    }
    if (zonal)
    {
        const volVectorField gRec(fvc::reconstruct(hf*mesh.magSf()));
        const boolList& zc = splitPtr->inZone();
        vectorField& gcf = gC.primitiveFieldRef();
        const scalarField& rhoc = rho.primitiveField();
        forAll(gcf, celli)
        {
            if (!zc[celli]) gcf[celli] = gRec[celli]/rhoc[celli];
        }
    }
    gC.correctBoundaryConditions();

    // gC is NOT applied here: it enters at the next step, through the
    // viscous shift and the permanent correction above

    p == p_rgh + rho*gh;
    if (p_rgh.needReference())
    {
        p += dimensionedScalar
        (
            "p",
            p.dimensions(),
            pRefValue - getRefCellValue(p, pRefCell)
        );
        p_rgh = p - rho*gh;
    }
}
