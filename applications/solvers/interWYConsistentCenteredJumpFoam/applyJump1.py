#!/usr/bin/env python3
"""Step 1 of interWYConsistentCenteredJumpFoam.

One-sided face density at the faces the interface crosses,

    1/rho_f = lambda/rho_s(P) + (1 - lambda)/rho_s(N),
    lambda  = psi_P/(psi_P - psi_N),

with psi the reconstructed distance function and s(.) the phase of a cell
by the sign of psi.  Away from those faces the harmonic mixture average is
kept, so the switch off reproduces the parent solver bit for bit.

Regenerates every touched file from its .pristine copy, so the edits never
stack.  Every replacement is asserted.
"""

import os
import shutil

D = ("/work1/u4485952/OpenFOAM/u4485952-v2506/geometricVoFCartesian/"
     "applications/solvers/interWYConsistentCenteredJumpFoam")
MAIN = "interWYConsistentCenteredJumpFoam.C"
TOUCHED = ["createFields.H", "pEqn.H", MAIN]


def restore(f):
    p = os.path.join(D, f)
    assert os.path.exists(p + ".pristine"), f"no pristine copy of {f}"
    shutil.copy2(p + ".pristine", p)
    return p


def sub(text, old, new, where, count=1):
    n = text.count(old)
    assert n == count, f"{n} occurrences of the anchor in {where}, expected {count}"
    return text.replace(old, new)


for f in TOUCHED:
    restore(f)

# ------------------------------------------------------------------ main .C
p = os.path.join(D, MAIN)
s = open(p).read()
s = sub(s, '#include "fvOptions.H"',
        '#include "fvOptions.H"\n'
        '#include "reconstructedDistanceFunction.H"\n'
        '#include "zoneDistribute.H"\n'
        '#include "syncTools.H"',
        "main .C includes")
open(p, "w").write(s)

# ------------------------------------------------------- createFields.H
p = restore("createFields.H")
s = open(p).read()
anchor = ("const bool staggerAlpha =\n"
          '    pimple.dict().getOrDefault<bool>("wyStaggerAlpha", true);')
s = sub(s, anchor, anchor + r'''

// Interfacial jump treatment of the face density in the pressure equation.
// At a face the interface crosses, the 1D-exact Poisson coefficient is the
// series resistance of the two phase segments,
//
//     1/rho_f = lambda/rho_s(P) + (1 - lambda)/rho_s(N),
//
// with lambda the interface position along the owner-neighbour segment and
// s(.) the phase of a cell.  The harmonic mixture average that the parent
// solver uses is this expression with lambda pinned to 1/2 and the mixture
// cell densities substituted for the phase densities
// (Vukcevic et al. 2017, Eqs. 39-49).  Default false = parent path.
const bool jumpFaceDensity
(
    pimple.dict().getOrDefault<Switch>("jumpFaceDensity", false)
);

// re-reconstruct the interface after the advection, so the distance
// function and the density field describe the same interface
const bool jumpReconstruct
(
    pimple.dict().getOrDefault<Switch>("jumpReconstruct", true)
);

autoPtr<reconstructedDistanceFunction> rdfPtr;
if (jumpFaceDensity)
{
    rdfPtr.reset(new reconstructedDistanceFunction(mesh));
    Info<< "jump: one-sided face density at interface faces"
        << " (jumpFaceDensity true, jumpReconstruct "
        << Switch(jumpReconstruct) << ")" << nl;
}
else
{
    Info<< "jump: face density is the harmonic mixture average"
        << " (jumpFaceDensity false, parent path)" << nl;
}''', "createFields.H switches")
open(p, "w").write(s)

# --------------------------------------------------------------- pEqn.H
p = restore("pEqn.H")
s = open(p).read()
s = sub(s,
        "    const surfaceScalarField rhof(1.0/fvc::interpolate(1.0/rho));",
        "    surfaceScalarField rhof(1.0/fvc::interpolate(1.0/rho));\n"
        "    if (jumpFaceDensity)\n"
        "    {\n"
        '        #include "jumpFaceDensity.H"\n'
        "    }",
        "pEqn.H face density")
open(p, "w").write(s)

# ----------------------------------------------------- jumpFaceDensity.H
open(os.path.join(D, "jumpFaceDensity.H"), "w").write(r'''// One-sided face density at the faces the interface crosses.
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

    auto oneSided = [r1, r2](const scalar psiP, const scalar psiN)
    {
        const scalar den = psiP - psiN;
        const scalar lambda = min(max(psiP/den, scalar(0)), scalar(1));
        const scalar rP = (psiP > 0) ? r1 : r2;
        const scalar rN = (psiN > 0) ? r1 : r2;
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
''')

# the report controls live with the other switches
p = os.path.join(D, "createFields.H")
s = open(p).read()
s = sub(s, "autoPtr<reconstructedDistanceFunction> rdfPtr;",
        "const bool jumpReport\n"
        "(\n"
        '    pimple.dict().getOrDefault<Switch>("jumpReport", true)\n'
        ");\n"
        "const label jumpReportInterval\n"
        "(\n"
        '    pimple.dict().getOrDefault<label>("jumpReportInterval", 100)\n'
        ");\n\n"
        "autoPtr<reconstructedDistanceFunction> rdfPtr;",
        "createFields.H report controls")
open(p, "w").write(s)

print("step 1 applied from pristine, all replacements asserted")
