/*---------------------------------------------------------------------------*\
    Part of interWYCenteredFoam (see the main file for the full
    description and license). Extracted for the conventional
    fragmented solver layout.
\*---------------------------------------------------------------------------*/

    // --- per-direction face-to-cell average machinery -----------------
    // (the approximate projection's cell correction; uniform Cartesian)
    const labelUList& own = mesh.owner();
    const labelUList& nei = mesh.neighbour();
    const surfaceVectorField& Sf = mesh.Sf();
    const surfaceScalarField& magSf = mesh.magSf();

    auto faceToCellAvg = [&]
    (
        const surfaceScalarField& gf,   // signed along +axis [m/s2]
        volVectorField& gcell
    )
    {
        gcell == dimensionedVector(dimAcceleration, Zero);
        vectorField& gc = gcell.primitiveFieldRef();

        forAll(own, facei)
        {
            const vector nHat = Sf[facei]/magSf[facei];
            for (direction d = 0; d < 3; ++d)
            {
                if (mag(nHat[d]) < 0.99) continue;
                const scalar v = gf[facei]*sign(nHat[d]);
                gc[own[facei]][d] += 0.5*v;
                gc[nei[facei]][d] += 0.5*v;
            }
        }
        forAll(mesh.boundary(), patchi)
        {
            const fvPatch& patch = mesh.boundary()[patchi];
            if (isA<emptyFvPatch>(patch)) continue;
            const labelUList& fc = patch.faceCells();
            const vectorField& pSf = Sf.boundaryField()[patchi];
            const scalarField& pMagSf = magSf.boundaryField()[patchi];
            const scalarField& pGf = gf.boundaryField()[patchi];
            forAll(fc, fi)
            {
                const vector nHat = pSf[fi]/pMagSf[fi];
                for (direction d = 0; d < 3; ++d)
                {
                    if (mag(nHat[d]) < 0.99) continue;
                    gc[fc[fi]][d] += 0.5*pGf[fi]*sign(nHat[d]);
                }
            }
        }
        gcell.correctBoundaryConditions();
    };

    // --- config echo --------------------------------------------------
    const bool staggerAlpha =
        pimple.dict().getOrDefault<bool>("wyStaggerAlpha", true);
    bool staggerDone = (runTime.value() > SMALL);
    scalar dtPrev = -1.0;

    Info<< "interWYCenteredFoam config:"
        << " reconstruction "
        << mesh.solverDict(alpha1.name())
              .getOrDefault<word>("reconstructionScheme", "plicRDF")
        << " staggerAlpha=" << staggerAlpha
        << " (single projection, averaged cell correction,"
        << " 1/interp(rho) face weights)"
        << nl << endl;

