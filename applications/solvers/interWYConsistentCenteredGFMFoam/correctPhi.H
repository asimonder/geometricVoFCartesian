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

    surfaceScalarField rhofCorr(1.0/fvc::interpolate(1.0/rho));
    if (gfmPressureJump)
    {
        // pairing rule: the topology-change projection must see the
        // same one-sided interface coefficient as pEqn, or every AMR
        // event reintroduces the smeared weights. pcorr is a jump-free
        // correction potential: coefficient only, no source.
        surfaceScalarField& gfmRhof = rhofCorr;
        const char* gfmCtx = "correctPhi";
        #include "gfmPressureJump.H"
    }
    const surfaceScalarField rAUfCorr("rAUf", dtDim/rhofCorr);

    if (pcorr.needReference())
    {
        fvc::makeRelative(phi, U);
        adjustPhi(phi, U, pcorr);
        fvc::makeAbsolute(phi, U);
    }

    mesh.setFluxRequired(pcorr.name());

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pcorrEqn
        (
            fvm::laplacian(rAUfCorr, pcorr) == fvc::div(phi)
        );

        pcorrEqn.setReference(pRefCell, 0);

        pcorrEqn.solve(pcorr.select(pimple.finalNonOrthogonalIter()));

        if (pimple.finalNonOrthogonalIter())
        {
            phi -= pcorrEqn.flux();
        }
    }
}
