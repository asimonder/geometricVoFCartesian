// 1. Setup Layer 0 (Interface)
volScalarField interfaceCells
(
    IOobject("interfaceCells", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
    mesh,
    dimensionedScalar(dimless, 0.0)
);

forAll(alpha1, celli)
{
    if (alpha1[celli] > interfaceThreshold && alpha1[celli] < (1.0 - interfaceThreshold))
    {
        interfaceCells[celli] = 1.0;
    }
}

// 2. EXPANSION LOOP (Topology Based)
if (nBufferLayers > 0)
{
    // Build/Get cell-to-cell connectivity (Standard OF Topology)
    const labelListList& cellCells = mesh.cellCells();

    for (label layer = 0; layer < nBufferLayers; ++layer)
    {
        // Snapshot: Read from 'old', write to 'interfaceCells'
        volScalarField oldLayer = interfaceCells;
        
        forAll(oldLayer, celli)
        {
            // If I am active, mark all my neighbors
            if (oldLayer[celli] > 0.5)
            {
                const labelList& neighbors = cellCells[celli];
                forAll(neighbors, ni)
                {
                    interfaceCells[neighbors[ni]] = 1.0;
                }
            }
        }
        
        // Sync processors (Ghost cells propagate the expansion)
        interfaceCells.correctBoundaryConditions();
    }
}

// 3. Create Flux Mask
surfaceScalarField interfaceMask
(
    IOobject("interfaceMask", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
    pos(fvc::interpolate(interfaceCells) - 0.001)
);

// 4. Split Fluxes
surfaceScalarField rhoPhiInter
(
    IOobject("rhoPhiInter", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
    rhoPhi * interfaceMask
);

surfaceScalarField rhoPhiBulk
(
    IOobject("rhoPhiBulk", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
    rhoPhi - rhoPhiInter
);
