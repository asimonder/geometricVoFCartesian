{
    // HARMONIC face density: the weight family of the droplet-passing
    // collocated solvers
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
        // cells
        af += (1.0 - fvc::interpolate(alpha1))
             *((airForce & mesh.Sf())/mesh.magSf())/rhof;
    }

    // Dev switch prognosticFlux (default off): propose the new flux
    // as the previous projected flux plus the step increment of the
    // cell field, so the checkerboard mismatch of the approximate cell
    // correction is never re-imported; full rebuild while the mesh
    // changes, where correctPhi has just restored phi
    const bool prognostic =
        pimple.dict().getOrDefault<bool>("prognosticFlux", false)
     && !mesh.changing();

    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        prognostic
      ? (phi
          + ((fvc::interpolate(U) - fvc::interpolate(Upre)) & mesh.Sf())
          + dtDim*af*mesh.magSf())
      : ((fvc::interpolate(U) & mesh.Sf()) + dtDim*af*mesh.magSf())
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

    // approximate-projection cell correction: the weight-free face
    // force density gathered to cells, divided by the CELL density.
    // reconstruct: general meshes (reduces to the two-face average on
    // uniform Cartesian grids); average: explicit axis average
    // (uniform axis-aligned meshes only).
    const surfaceScalarField gf(af - (rAUf/dtDim)*fvc::snGrad(p_rgh));
    const surfaceScalarField hf(rhof*gf);
    if (cellCorrection == "average")
    {
        faceToCellAvg(hf, gC);
        vectorField& gcf = gC.primitiveFieldRef();
        const scalarField& rhoc = rho.primitiveField();
        forAll(gcf, celli)
        {
            gcf[celli] /= rhoc[celli];
        }
    }
    else
    {
        const volVectorField gRec(fvc::reconstruct(hf*mesh.magSf()));
        vectorField& gcf = gC.primitiveFieldRef();
        const scalarField& rhoc = rho.primitiveField();
        forAll(gcf, celli)
        {
            gcf[celli] = gRec[celli]/rhoc[celli];
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
