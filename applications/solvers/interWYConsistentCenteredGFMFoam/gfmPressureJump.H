// Phase-1 Ghost Fluid Method, coefficient-only form: the one-sided
// subface face density at interface faces (Vukcevic et al. 2017,
// Eqs. 39-42 marking/position, Eqs. 46-47/52-53 coefficient). No jump
// source and no dynamic-pressure redefinition: gravity remains in the
// solver's balanced-force form, so its discrete work stays tied to the
// geometric alpha transport (the N0128 energy-injection lesson,
// 2026-07-25).
//
// Included from pEqn.H and correctPhi.H so both Laplacians share the
// one-sided interface weights.
//
// Expects in scope:
//   surfaceScalarField& gfmRhof    face density; interface faces are
//                                  overwritten with the subface
//                                  harmonic of the PURE densities,
//                                  1/rho_f = lam/rhoP,pure
//                                          + (1-lam)/rhoN,pure.
//                                  Light-side-scaled like the
//                                  production harmonic, so the
//                                  correction pairing and its
//                                  stability bound survive; sharp in
//                                  the sense that the smeared cell
//                                  values are replaced by pure
//                                  densities at the geometric
//                                  interface position. NOT the
//                                  Vukcevic Eq. 46-47 coefficient:
//                                  that one is water-scaled and is
//                                  measured incompatible with this
//                                  cycle's correction (2026-07-25,
//                                  injection with acceleration
//                                  delivery, blowup with force
//                                  delivery).
//   const char* gfmCtx             log-line context tag
//
// Classification: wet = alpha1 > 0.5. A cell at exactly 0.5 marks no
// interface face; the Eq. 39 product criterion is strict.
{
    // lambda must come from the interface the projection acts on:
    // re-reconstruct against the current alpha (no-op if the curvature
    // machinery already did this step)
    advector.surf().reconstruct(false);

    const volScalarField& gfmRDF =
        mesh.lookupObject<volScalarField>("RDF");
    const volScalarField& gfmCDL =
        mesh.lookupObject<volScalarField>("cellDistLevel");

    // stale-halo guard (the reconstruction-normal parallel lesson):
    // sync coupled values before any patchNeighbourField read
    const_cast<volScalarField&>(gfmRDF).correctBoundaryConditions();
    const_cast<volScalarField&>(gfmCDL).correctBoundaryConditions();

    const scalar gfmRhoW = rho1.value();
    const scalar gfmRhoA = rho2.value();

    label gfmNIn = 0;        // interface faces, internal
    label gfmNCp = 0;        // interface face-sides on coupled patches
    label gfmNClamp = 0;     // lambda clamped into [0,1]
    label gfmNFall = 0;      // alpha-ratio fallbacks (Eq. 41)
    label gfmNPhys = 0;      // interface contact on physical patches

    // one face: subface position -> one-sided coefficient
    auto gfmFace = [&]
    (
        const scalar aP, const scalar aN,
        const scalar rdfP, const scalar rdfN,
        const scalar cdlP, const scalar cdlN
    ) -> scalar
    {
        scalar lam;
        const scalar den = rdfP - rdfN;
        if
        (
            cdlP >= 0 && cdlN >= 0
         && mag(den) > max(VSMALL, 1e-8*(mag(rdfP) + mag(rdfN)))
        )
        {
            lam = rdfP/den;                      // Eq. 40, RDF variant
        }
        else
        {
            lam = (aP - 0.5)/(aP - aN);          // Eq. 41 fallback
            ++gfmNFall;
        }
        if (lam < 0.0) { lam = 0.0; ++gfmNClamp; }
        else if (lam > 1.0) { lam = 1.0; ++gfmNClamp; }

        const bool wetP = (aP > 0.5);
        const scalar rhoP = wetP ? gfmRhoW : gfmRhoA;
        const scalar rhoN = wetP ? gfmRhoA : gfmRhoW;
        return 1.0/(lam/rhoP + (1.0 - lam)/rhoN);
    };

    {
        const labelUList& own = mesh.owner();
        const labelUList& nei = mesh.neighbour();
        const scalarField& aIn = alpha1.primitiveField();
        const scalarField& rdfIn = gfmRDF.primitiveField();
        const scalarField& cdlIn = gfmCDL.primitiveField();

        scalarField& rhofIn = gfmRhof.primitiveFieldRef();

        forAll(own, facei)
        {
            const scalar aP = aIn[own[facei]];
            const scalar aN = aIn[nei[facei]];
            if ((aP - 0.5)*(aN - 0.5) >= 0.0) continue;   // Eq. 39
            ++gfmNIn;

            rhofIn[facei] = gfmFace
            (
                aP, aN,
                rdfIn[own[facei]], rdfIn[nei[facei]],
                cdlIn[own[facei]], cdlIn[nei[facei]]
            );
        }
    }

    forAll(mesh.boundary(), patchi)
    {
        const fvPatch& patch = mesh.boundary()[patchi];
        if (isA<emptyFvPatch>(patch)) continue;

        const fvPatchScalarField& aBf = alpha1.boundaryField()[patchi];

        if (!patch.coupled())
        {
            // no interface treatment on physical patches in phase 1;
            // count contact so nothing passes silently
            const scalarField aPc(aBf.patchInternalField());
            forAll(aPc, fi)
            {
                if ((aPc[fi] - 0.5)*(aBf[fi] - 0.5) < 0.0) ++gfmNPhys;
            }
            continue;
        }

        const scalarField aPc(aBf.patchInternalField());
        const scalarField aNc(aBf.patchNeighbourField());
        const scalarField rdfPc
        (
            gfmRDF.boundaryField()[patchi].patchInternalField()
        );
        const scalarField rdfNc
        (
            gfmRDF.boundaryField()[patchi].patchNeighbourField()
        );
        const scalarField cdlPc
        (
            gfmCDL.boundaryField()[patchi].patchInternalField()
        );
        const scalarField cdlNc
        (
            gfmCDL.boundaryField()[patchi].patchNeighbourField()
        );

        scalarField& rhofP = gfmRhof.boundaryFieldRef()[patchi];

        forAll(aPc, fi)
        {
            if ((aPc[fi] - 0.5)*(aNc[fi] - 0.5) >= 0.0) continue;
            ++gfmNCp;

            rhofP[fi] = gfmFace
            (
                aPc[fi], aNc[fi], rdfPc[fi], rdfNc[fi],
                cdlPc[fi], cdlNc[fi]
            );
        }
    }

    reduce(gfmNIn, sumOp<label>());
    reduce(gfmNCp, sumOp<label>());
    reduce(gfmNClamp, sumOp<label>());
    reduce(gfmNFall, sumOp<label>());
    reduce(gfmNPhys, sumOp<label>());

    // coupled face-sides are seen twice (cyclic halves, processor twins)
    Info<< "GFM(" << gfmCtx << "): "
        << gfmNIn + gfmNCp/2 << " interface faces, "
        << gfmNClamp << " lambda clamps, "
        << gfmNFall << " alpha fallbacks";
    if (gfmNPhys > 0)
    {
        Info<< ", " << gfmNPhys
            << " interface contacts on physical patches (no"
               " treatment)";
    }
    Info<< endl;
}
