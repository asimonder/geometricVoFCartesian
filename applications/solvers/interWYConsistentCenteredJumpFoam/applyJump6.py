#!/usr/bin/env python3
"""Step 6 of interWYConsistentCenteredJumpFoam.

One-sided cell correction next to the interface.

For each direction independently, a cell normally averages the two faces
whose normal lies along that direction.  Where one of those two faces
touches a mixed cell it carries the band cancellation residual, which the
air cell then divides by its own small density.  This step drops that
face from the stencil and takes the acceleration from the remaining one,
rescaled so the weight is unchanged.  Where both faces are clean, or both
touch a mixed cell, the parent average is kept exactly.

No force is changed.  Only the set of faces a cell listens to changes, and
the test is alpha in the neighbouring cell, so no second interface
descriptor enters.

At rest gf vanishes on every face, so the hydrostatic balance cannot see
this change.

Runs applyJump5.py first, so all steps are regenerated from pristine.
"""

import os
import runpy

D = ("/work1/u4485952/OpenFOAM/u4485952-v2506/geometricVoFCartesian/"
     "applications/solvers/interWYConsistentCenteredJumpFoam")

runpy.run_path(os.path.join(D, "applyJump5.py"))


def sub(text, old, new, where, count=1):
    n = text.count(old)
    assert n == count, f"{n} occurrences of the anchor in {where}, expected {count}"
    return text.replace(old, new)


# ------------------------------------------------------- createFields.H
p = os.path.join(D, "createFields.H")
s = open(p).read()

anchor = "const bool staggerAlpha ="
s = sub(s, anchor, r'''// One-sided cell correction beside the interface.  Per direction, a face
// touching a mixed cell is dropped from the cell's stencil and the
// remaining face carries the whole weight.  Identical to the parent
// wherever a direction's two faces are both clean or both touch a mixed
// cell.  Default false = parent path.
const bool cellCorrectionOneSided
(
    pimple.dict().getOrDefault<Switch>("cellCorrectionOneSided", false)
);
const scalar cellCorrectionAlphaTol
(
    pimple.dict().getOrDefault<scalar>("cellCorrectionAlphaTol", 1e-6)
);
if (cellCorrectionOneSided)
{
    Info<< "jump: one-sided cell correction beside the interface"
        << " (cellCorrectionOneSided true, alphaTol "
        << cellCorrectionAlphaTol << ")" << nl;
}

// per-direction face-to-cell average that drops the faces touching a
// mixed cell.  Same 0.5 weights as the parent, rescaled by the number of
// faces actually used, so a direction with two clean faces is identical.
auto faceToCellAvgOneSided =
[&mesh, &alpha1, &cellCorrectionAlphaTol]
(
    const surfaceScalarField& gfld,
    volVectorField& gcell
)
{
    const labelUList& own = mesh.owner();
    const labelUList& nei = mesh.neighbour();
    const surfaceVectorField& Sf = mesh.Sf();
    const surfaceScalarField& magSf = mesh.magSf();
    const scalarField& a = alpha1.primitiveField();
    const scalar aTol = cellCorrectionAlphaTol;

    gcell == dimensionedVector(dimAcceleration, Zero);
    vectorField& gc = gcell.primitiveFieldRef();
    vectorField cleanSum(gc.size(), Zero);
    vectorField nAll(gc.size(), Zero);
    vectorField nClean(gc.size(), Zero);

    auto isMixed = [&a, aTol](const label c)
    {
        return a[c] > aTol && a[c] < 1 - aTol;
    };
    auto add =
    [&gc, &cleanSum, &nAll, &nClean]
    (const label c, const direction d, const scalar v, const bool clean)
    {
        gc[c][d] += 0.5*v;
        nAll[c][d] += 1;
        if (clean)
        {
            cleanSum[c][d] += 0.5*v;
            nClean[c][d] += 1;
        }
    };

    forAll(own, facei)
    {
        const vector nHat = Sf[facei]/magSf[facei];
        const bool clean = !isMixed(own[facei]) && !isMixed(nei[facei]);
        for (direction d = 0; d < 3; ++d)
        {
            if (mag(nHat[d]) < 0.99) continue;
            const scalar v = gfld[facei]*sign(nHat[d]);
            add(own[facei], d, v, clean);
            add(nei[facei], d, v, clean);
        }
    }

    forAll(mesh.boundary(), patchi)
    {
        const fvPatch& patch = mesh.boundary()[patchi];
        if (isA<emptyFvPatch>(patch)) continue;

        const labelUList& fc = patch.faceCells();
        const vectorField& pSf = Sf.boundaryField()[patchi];
        const scalarField& pMagSf = magSf.boundaryField()[patchi];
        const scalarField& pGf = gfld.boundaryField()[patchi];

        tmp<scalarField> taNei;
        if (patch.coupled())
        {
            taNei = alpha1.boundaryField()[patchi].patchNeighbourField();
        }

        forAll(fc, fi)
        {
            bool clean = !isMixed(fc[fi]);
            if (clean && patch.coupled())
            {
                const scalar an = taNei()[fi];
                clean = !(an > aTol && an < 1 - aTol);
            }
            const vector nHat = pSf[fi]/pMagSf[fi];
            for (direction d = 0; d < 3; ++d)
            {
                if (mag(nHat[d]) < 0.99) continue;
                add(fc[fi], d, pGf[fi]*sign(nHat[d]), clean);
            }
        }
    }

    forAll(gc, celli)
    {
        for (direction d = 0; d < 3; ++d)
        {
            const scalar na = nAll[celli][d];
            const scalar nc = nClean[celli][d];
            if (na > 0 && nc > 0 && nc < na)
            {
                gc[celli][d] = cleanSum[celli][d]*(na/nc);
            }
        }
    }
    gcell.correctBoundaryConditions();
};

''' + anchor, "createFields.H step 6 lambda")
open(p, "w").write(s)

# --------------------------------------------------------------- pEqn.H
p = os.path.join(D, "pEqn.H")
s = open(p).read()
s = sub(s,
        "        faceToCellAvg(hf, gC);\n"
        "        vectorField& gcf = gC.primitiveFieldRef();",
        "        if (cellCorrectionOneSided)\n"
        "        {\n"
        "            faceToCellAvgOneSided(hf, gC);\n"
        "        }\n"
        "        else\n"
        "        {\n"
        "            faceToCellAvg(hf, gC);\n"
        "        }\n"
        "        vectorField& gcf = gC.primitiveFieldRef();",
        "pEqn.H one-sided branch")
open(p, "w").write(s)

print("step 6 applied from pristine, all replacements asserted")
