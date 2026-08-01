// One-sided face density at the faces the interface crosses.
//
//     1/rho_f = lambda/rho_s(P) + (1 - lambda)/rho_s(N)
//     lambda  = psi_P/(psi_P - psi_N)
//
// psi is the reconstructed distance function, positive in phase 1.  Faces
// with no sign change keep the harmonic mixture average built above.
// Coupled patches are treated exactly like internal faces, so a cyclic or
// processor boundary has no effect on the result.
{
    reconstructionSchemes& surf = advector.surf();
    if (jumpReconstruct)
    {
        surf.reconstruct();
    }

    const volVectorField& iCentre = surf.centre();
    const volVectorField& iNormal = surf.normal();

    boolList interfaceCells(mesh.nCells(), false);
    forAll(iNormal, celli)
    {
        interfaceCells[celli] = (mag(iNormal[celli]) != 0);
    }

    // a two-layer band covers both cells of every crossed face
    rdfPtr->markCellsNearSurf(interfaceCells, 2);
    const boolList& band = rdfPtr->nextToInterface();

    // stale values outside the band would survive from the previous step
    volScalarField& psiRef = *rdfPtr;
    psiRef = dimensionedScalar(dimLength, Zero);

    zoneDistribute& exchange = zoneDistribute::New(mesh);
    const volScalarField& psi =
        rdfPtr->constructRDF(band, iCentre, iNormal, exchange, true);

    scalarField bandCell(mesh.nCells(), Zero);
    forAll(band, celli)
    {
        bandCell[celli] = band[celli] ? 1.0 : 0.0;
    }
    scalarField psiNei, bandNei;
    syncTools::swapBoundaryCellList(mesh, psi.primitiveField(), psiNei);
    syncTools::swapBoundaryCellList(mesh, bandCell, bandNei);

    const scalar r1 = rho1.value();
    const scalar r2 = rho2.value();

    const bool qv = (jumpFaceDensityModel == "queutey");
    auto oneSided = [r1, r2, qv](const scalar psiP, const scalar psiN)
    {
        const scalar rP = (psiP > 0) ? r1 : r2;
        const scalar rN = (psiN > 0) ? r1 : r2;
        if (qv)
        {
            // Queutey and Visonneau (2007) Eq. (36b) on a uniform
            // orthogonal mesh with pure phases either side
            return 0.5*(rP + rN);
        }
        const scalar den = psiP - psiN;
        const scalar lambda = min(max(psiP/den, scalar(0)), scalar(1));
        return scalar(1)/(lambda/rP + (scalar(1) - lambda)/rN);
    };

    label nJump = 0;
    {
        const labelUList& own = mesh.owner();
        const labelUList& nei = mesh.neighbour();
        scalarField& rhofi = rhof.primitiveFieldRef();

        forAll(own, facei)
        {
            const label P = own[facei];
            const label N = nei[facei];
            if (!band[P] || !band[N]) continue;

            const scalar psiP = psi[P];
            const scalar psiN = psi[N];
            if (psiP*psiN > 0) continue;
            if (mag(psiP - psiN) < SMALL) continue;

            rhofi[facei] = oneSided(psiP, psiN);
            ++nJump;
        }
    }

    forAll(mesh.boundary(), patchi)
    {
        const fvPatch& pp = mesh.boundary()[patchi];
        if (!pp.coupled()) continue;

        const labelUList& fc = pp.faceCells();
        const label start = pp.patch().start() - mesh.nInternalFaces();
        scalarField& rhofp = rhof.boundaryFieldRef()[patchi];

        forAll(fc, fi)
        {
            const label P = fc[fi];
            const scalar psiP = psi[P];
            const scalar psiN = psiNei[start + fi];
            if (!band[P] || bandNei[start + fi] < 0.5) continue;
            if (psiP*psiN > 0) continue;
            if (mag(psiP - psiN) < SMALL) continue;

            rhofp[fi] = oneSided(psiP, psiN);
            ++nJump;
        }
    }

    if (jumpReport && (runTime.timeIndex() % jumpReportInterval == 0))
    {
        Info<< "jump: one-sided density on "
            << returnReduce(nJump, sumOp<label>()) << " faces, rho_f in ["
            << gMin(rhof.primitiveField()) << ", "
            << gMax(rhof.primitiveField()) << "] kg/m3" << endl;
    }
}
