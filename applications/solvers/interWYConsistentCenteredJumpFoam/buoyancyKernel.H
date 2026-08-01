// Buoyancy distributed by a regularised Heaviside of the distance
// function.  Replaces snGrad(rho) in the gravity face force only; the
// transported density, the face density and the transport itself are
// untouched.
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

    rdfPtr->markCellsNearSurf(interfaceCells, 2);
    const boolList& band = rdfPtr->nextToInterface();

    volScalarField& psiRef = *rdfPtr;
    psiRef = dimensionedScalar(dimLength, Zero);

    zoneDistribute& exchange = zoneDistribute::New(mesh);
    const volScalarField& psi =
        rdfPtr->constructRDF(band, iCentre, iNormal, exchange, true);

    // one global cell size; the sweep machinery already requires a
    // uniform lattice where the interface lives
    const scalar hCell =
        1.0/gMax(mesh.surfaceInterpolation::deltaCoeffs().primitiveField());
    const scalar eps = max(buoyancyKernelWidth*hCell, SMALL);
    const scalar pi = Foam::constant::mathematical::pi;

    // cell-centred regularised Heaviside; outside the band the phase is
    // taken from alpha so that the telescoping closes at the band edge
    scalarField H(mesh.nCells());
    const scalarField& a = alpha1.primitiveField();
    forAll(H, celli)
    {
        if (!band[celli])
        {
            H[celli] = (a[celli] > 0.5) ? 1.0 : 0.0;
            continue;
        }
        const scalar r = psi[celli]/eps;
        if (r <= -1.0)      H[celli] = 0.0;
        else if (r >= 1.0)  H[celli] = 1.0;
        else                H[celli] = 0.5*(1.0 + r + Foam::sin(pi*r)/pi);
    }

    scalarField HNei;
    syncTools::swapBoundaryCellList(mesh, H, HNei);

    const scalar dRho = rho1.value() - rho2.value();
    surfaceScalarField& sg = tsnGradRho.ref();

    {
        const labelUList& own = mesh.owner();
        const labelUList& nei = mesh.neighbour();
        const surfaceScalarField& dc = mesh.deltaCoeffs();
        scalarField& sgi = sg.primitiveFieldRef();

        forAll(own, facei)
        {
            sgi[facei] = dRho*(H[nei[facei]] - H[own[facei]])*dc[facei];
        }
    }

    forAll(mesh.boundary(), patchi)
    {
        const fvPatch& pp = mesh.boundary()[patchi];
        if (!pp.coupled()) continue;

        const labelUList& fc = pp.faceCells();
        const label start = pp.patch().start() - mesh.nInternalFaces();
        const scalarField& dc = mesh.deltaCoeffs().boundaryField()[patchi];
        scalarField& sgp = sg.boundaryFieldRef()[patchi];

        forAll(fc, fi)
        {
            sgp[fi] = dRho*(HNei[start + fi] - H[fc[fi]])*dc[fi];
        }
    }

    if (jumpReport && (runTime.timeIndex() % jumpReportInterval == 0))
    {
        label nb = 0;
        forAll(band, celli) if (band[celli]) ++nb;
        Info<< "jump: buoyancy kernel eps " << eps << " m over "
            << returnReduce(nb, sumOp<label>()) << " band cells" << endl;
    }
}
