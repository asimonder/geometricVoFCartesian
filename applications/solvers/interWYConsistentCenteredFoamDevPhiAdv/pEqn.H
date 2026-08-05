{
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

    // DevPhiAdv face-carried momentum predictor: phi advances by its
    // own Lax-Wendroff face advection plus the face-interpolated
    // viscous increment; the cell velocity enters only through that
    // increment, severing the cell-noise feedback. Faces without a
    // complete lattice stencil keep the classical predictor.
    surfaceScalarField phiAdv
    (
        IOobject("phiAdv", runTime.timeName(), mesh,
                 IOobject::NO_READ, IOobject::NO_WRITE,
                 IOobject::NO_REGISTER),
        (fvc::interpolate(UpreVisc) & mesh.Sf())
    );
    {
        const scalar dtv = runTime.deltaTValue();
        const surfaceScalarField& magSfl = mesh.magSf();
        const scalarField& phif = phi.primitiveField();
        const scalarField& msf = magSfl.primitiveField();
        scalarField& pAdv = phiAdv.primitiveFieldRef();

        for (label fi = 0; fi < mesh.nInternalFaces(); ++fi)
        {
            if (!faceOkAdv[fi]) continue;
            const label d = faceAxis[fi];
            const scalar ud = phif[fi]/msf[fi];
            scalar A = 0;
            for (direction e = 0; e < 3; ++e)
            {
                if (axEmpty[e]) continue;
                const label fm = fAdj[fi][2*e];
                const label fp = fAdj[fi][2*e+1];
                const scalar um = phif[fm]/msf[fm];
                const scalar up = phif[fp]/msf[fp];
                scalar ue;
                if (label(e) == d)
                {
                    ue = ud;
                }
                else
                {
                    ue = 0.25*(
                        phif[uAdj[fi][4*e]]/msf[uAdj[fi][4*e]]
                      + phif[uAdj[fi][4*e+1]]/msf[uAdj[fi][4*e+1]]
                      + phif[uAdj[fi][4*e+2]]/msf[uAdj[fi][4*e+2]]
                      + phif[uAdj[fi][4*e+3]]/msf[uAdj[fi][4*e+3]]);
                }
                const scalar h = hLat[e];
                // Lax-Wendroff along axis e: centred slope plus the
                // stabilising second-difference correction
                A -= ue*(up - um)/(2.0*h);
                A += 0.5*dtv*ue*ue*(up - 2.0*ud + um)/(h*h);
            }
            pAdv[fi] = phif[fi] + dtv*A*msf[fi];
        }
    }

    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        phiAdv
      + (fvc::interpolate(U - UpreVisc) & mesh.Sf())
      + dtDim*af*mesh.magSf()
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
    // face force density, axis-averaged, divided by the CELL density
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
