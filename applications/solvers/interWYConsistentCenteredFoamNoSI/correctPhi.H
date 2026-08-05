// Flux restoration after a topology change, replacing the library
// CorrectPhi call. Two defects of the library route, both measured on
// the N0032 AMR reproducer:
//
//  1. CorrectPhi executes pcorrEqn.setReference(0, 0) unconditionally,
//     pinning the LOCAL cell 0 of every rank; each pinned cell retains
//     a divergence diag*pcorr that the solver residual cannot see.
//     Cure: the solver's single global reference cell.
//
//  2. The engine's recreated fluxes after unrefinement can carry a
//     nonzero GLOBAL divergence integral (inconsistent coupled-seam
//     recreation). No projection can remove a global imbalance; it
//     parks at the reference cell and converts into volume fraction
//     when that cell is interface-adjacent. Cure: rebuild the flux
//     from the boundary-synced velocity, which gives antisymmetric
//     seam fluxes by construction.
//
// The rebuild is CONDITIONAL: the projection preserves the global
// divergence integral, so compatibility is decided on the mapped flux
// before projecting. Compatible events (refinement-dominated; the
// large majority) keep the mapped pEqn flux structure and are only
// projected; an unconditional every-event rebuild perturbed the
// band-edge transport enough to trip the interface confinement guard
// at refineInterval 1 (N0256 devFix job 985467, t = 0.656).
{
    correctUphiBCs(U, phi);

    {
        const volScalarField divPhiMap(fvc::div(phi));
        const scalar imb =
            mag(gSum(divPhiMap.primitiveField()*mesh.V().field()));
        const scalar l1 =
            gSum(mag(divPhiMap.primitiveField())*mesh.V().field());

        if (imb > 1e-8*(l1 + VSMALL))
        {
            Info<< "wy: mapped flux global divergence imbalance " << imb
                << " (L1 " << l1 << "); rebuilding flux from velocity"
                << nl;
            U.correctBoundaryConditions();
            phi = fvc::interpolate(U) & mesh.Sf();
        }
    }

    wordList pcorrTypes
    (
        p_rgh.boundaryField().size(),
        fvPatchFieldBase::zeroGradientType()
    );

    forAll(p_rgh.boundaryField(), patchi)
    {
        if (p_rgh.boundaryField()[patchi].fixesValue())
        {
            pcorrTypes[patchi] = fixedValueFvPatchScalarField::typeName;
        }
    }

    volScalarField pcorr
    (
        IOobject
        (
            "pcorr",
            runTime.timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar(p_rgh.dimensions(), Zero),
        pcorrTypes
    );

    const surfaceScalarField rhofCorr(1.0/fvc::interpolate(1.0/rho));
    const surfaceScalarField rAUfCorr("rAUf", dtDim/rhofCorr);

    if (pcorr.needReference())
    {
        fvc::makeRelative(phi, U);
        adjustPhi(phi, U, pcorr);
        fvc::makeAbsolute(phi, U);
    }

    mesh.setFluxRequired(pcorr.name());

    // pcorr runs only at topology changes, so it carries its own
    // non-orthogonal corrector count: hanging-node faces created by the
    // event need more explicit correction passes than the per-step
    // p_rgh solve (measured ~14x residual reduction per pass), and the
    // cost here is negligible. Default = the PIMPLE count, unchanged
    // behaviour for existing cases.
    const label nNonOrthPcorr =
        pimple.dict().getOrDefault<label>
        (
            "nNonOrthogonalCorrectorsPcorr",
            pimple.nNonOrthCorr()
        );

    for (label nonOrth = 0; nonOrth <= nNonOrthPcorr; ++nonOrth)
    {
        fvScalarMatrix pcorrEqn
        (
            fvm::laplacian(rAUfCorr, pcorr) == fvc::div(phi)
        );

        pcorrEqn.setReference(pRefCell, 0);

        pcorrEqn.solve(pcorr.select(nonOrth == nNonOrthPcorr));

        if (nonOrth == nNonOrthPcorr)
        {
            phi -= pcorrEqn.flux();
        }
    }
}
