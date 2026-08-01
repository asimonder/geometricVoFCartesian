{
    // HARMONIC face density: the weight family of the droplet-passing
    // collocated solvers (rAUf = interp(dt/rho) in interIsoCartFoam);
    // GFM branch: interface faces carry the one-sided subface density
    // (Vukcevic et al. 2017, Eqs. 46-47, 52-53); gravity stays in the
    // balanced-force form below and shares the modified weights, so
    // its work remains tied to the alpha transport (no relocatable
    // jump source, the N0128 energy-injection lesson)
    surfaceScalarField rhof(1.0/fvc::interpolate(1.0/rho));
    if (gfmPressureJump)
    {
        surfaceScalarField& gfmRhof = rhof;
        const char* gfmCtx = "pEqn";
        #include "gfmPressureJump.H"
    }
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
        af += (scalar(1) - fvc::interpolate(alpha1))
             *((airForce & mesh.Sf())/mesh.magSf())/rhof;
    }

    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        (fvc::interpolate(U) & mesh.Sf()) + dtDim*af*mesh.magSf()
    );

    constrainPressure(p_rgh, U, phiHbyA, rAUf);

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix p_rghEqn
        (
            fvm::laplacian(rAUf, p_rgh) == fvc::div(phiHbyA)
        );
        p_rghEqn.setReference
        (
            pRefCell, getRefCellValue(p_rgh, pRefCell)
        );
        p_rghEqn.solve(p_rgh.select(pimple.finalNonOrthogonalIter()));
        if (pimple.finalNonOrthogonalIter())
        {
            phi = phiHbyA - p_rghEqn.flux();
        }
    }
    #include "continuityErrs.H"

    // averaged cell correction (approximate projection): the weight-free
    // face force density, axis-averaged, divided by the CELL density.
    // GFM branch: marked faces are excluded here and deliver the
    // continuous face acceleration instead (Vukcevic et al. 2017,
    // Eq. 44 discretised as the interface-corrected cell gradient,
    // Eqs. 62-65); a force/rho_cell delivery at a jump face amplifies
    // perturbations by rho_eff/rho_cell (measured ~100/step, N32 flat)
    const surfaceScalarField gf(af - (rAUf/dtDim)*fvc::snGrad(p_rgh));
    const surfaceScalarField hf(rhof*gf);
    if (gfmAirStencil)
    {
        faceToCellAvgOneSided(hf, alpha1, gC);
    }
    else
    {
        faceToCellAvg(hf, gC);
    }
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

    U += dtDim*gC;
    U.correctBoundaryConditions();

    fvOptions.correct(U);

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
