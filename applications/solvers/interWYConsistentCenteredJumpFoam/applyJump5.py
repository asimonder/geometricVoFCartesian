#!/usr/bin/env python3
"""Step 5 of interWYConsistentCenteredJumpFoam.

Buoyancy distributed by a regularised Heaviside of the distance function
instead of by the volume-fraction differences.

    H(r) = 0, r <= -1;  1, r >= 1;  (1 + r + sin(pi r)/pi)/2 otherwise
    snGrad(rho)_f  ->  (rho1 - rho2) (H_N - H_P) / |d_f|

with r = psi/eps, psi the reconstructed distance function and eps the
kernel width.  H is the smeared Heaviside of Sussman, Smereka and Osher
(1994), the standard regularisation whose derivative is the smeared delta
built to make an interfacial force independent of the sub-cell interface
position.

Conservation.  H_N - H_P telescopes along any line to H(top) - H(bottom),
which is 1 across the interface, so the total buoyancy per column is
exactly (rho1 - rho2) as before.  H is symmetric about psi = 0, so the
first moment stays at the interface.  Outside the band H is taken from the
sign of alpha - 1/2, which is 0 or 1 there, so the telescoping closes.

Mass transport is untouched.  Only the momentum source changes, phi is
still exactly projected, and the sweeps conserve alpha for any solenoidal
flux.

Runs applyJump4.py first, so all steps are regenerated from pristine.
"""

import os
import runpy

D = ("/work1/u4485952/OpenFOAM/u4485952-v2506/geometricVoFCartesian/"
     "applications/solvers/interWYConsistentCenteredJumpFoam")

runpy.run_path(os.path.join(D, "applyJump4.py"))


def sub(text, old, new, where, count=1):
    n = text.count(old)
    assert n == count, f"{n} occurrences of the anchor in {where}, expected {count}"
    return text.replace(old, new)


# ------------------------------------------------------- createFields.H
p = os.path.join(D, "createFields.H")
s = open(p).read()
anchor = "autoPtr<reconstructedDistanceFunction> rdfPtr;"
s = sub(s, anchor, r'''// Buoyancy distributed by a regularised Heaviside of the distance
// function rather than by the volume-fraction differences.  The total and
// the first moment per column are unchanged; what changes is that the
// distribution now varies continuously as the interface crosses a cell
// instead of jumping.  Sussman, Smereka and Osher (1994) smeared
// Heaviside.  Default false = parent path.
const bool buoyancyKernel
(
    pimple.dict().getOrDefault<Switch>("buoyancyKernel", false)
);
const scalar buoyancyKernelWidth
(
    pimple.dict().getOrDefault<scalar>("buoyancyKernelWidth", 1.0)
);
if (buoyancyKernel)
{
    Info<< "jump: buoyancy distributed by a regularised Heaviside, width "
        << buoyancyKernelWidth << " cells" << nl;
}

''' + anchor, "createFields.H step 5 switch")
s = sub(s, "if (jumpFaceDensity)\n{\n    rdfPtr.reset(new reconstructedDistanceFunction(mesh));",
        "if (jumpFaceDensity || buoyancyKernel)\n{\n    rdfPtr.reset(new reconstructedDistanceFunction(mesh));",
        "createFields.H rdf construction")
open(p, "w").write(s)

# --------------------------------------------------------------- pEqn.H
p = os.path.join(D, "pEqn.H")
s = open(p).read()
s = sub(s,
        "    surfaceScalarField af(-ghf*fvc::snGrad(rho)/rhof);",
        "    tmp<surfaceScalarField> tsnGradRho(fvc::snGrad(rho));\n"
        "    if (buoyancyKernel)\n"
        "    {\n"
        '        #include "buoyancyKernel.H"\n'
        "    }\n"
        "    surfaceScalarField af(-ghf*tsnGradRho()/rhof);",
        "pEqn.H buoyancy hook")
open(p, "w").write(s)

# ----------------------------------------------------- buoyancyKernel.H
open(os.path.join(D, "buoyancyKernel.H"), "w").write(r'''// Buoyancy distributed by a regularised Heaviside of the distance
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
''')

print("step 5 applied from pristine, all replacements asserted")
