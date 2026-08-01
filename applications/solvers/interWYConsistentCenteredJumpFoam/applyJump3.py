#!/usr/bin/env python3
"""Step 3 of interWYConsistentCenteredJumpFoam.

Cell correction from the face acceleration instead of the face force.

The parent averages the face force density hf = rho_f gf over the two
opposite faces of a cell and divides by the mixture cell density,

    gC = avg(hf)/rho_cell .

In a mixed cell the two faces carry forces that differ by the density
ratio, and the mixture density belongs to neither phase, so the resulting
acceleration belongs to neither phase either.  The face acceleration
gf = af - snGrad(p_rgh)/rho_f is the material acceleration, which is
continuous across the interface, so averaging it is well posed:

    gC = avg(gf) .

Where the density is constant rho_f = rho_cell and the two forms are
identical, so the change is confined to the interface band.  At rest
gf = 0 on every face, so the hydrostatic balance is untouched.

Runs applyJump2.py first, so all steps are regenerated from pristine.
"""

import os
import runpy

D = ("/work1/u4485952/OpenFOAM/u4485952-v2506/geometricVoFCartesian/"
     "applications/solvers/interWYConsistentCenteredJumpFoam")

runpy.run_path(os.path.join(D, "applyJump2.py"))


def sub(text, old, new, where, count=1):
    n = text.count(old)
    assert n == count, f"{n} occurrences of the anchor in {where}, expected {count}"
    return text.replace(old, new)


# ------------------------------------------------------- createFields.H
p = os.path.join(D, "createFields.H")
s = open(p).read()
anchor = "autoPtr<reconstructedDistanceFunction> rdfPtr;"
s = sub(s, anchor, r'''// Cell correction from the face acceleration rather than the face force.
// The parent averages hf = rho_f gf and divides by the mixture cell
// density.  In a mixed cell that mixes two face forces separated by the
// density ratio and divides by a density belonging to neither phase.
// The face acceleration gf is the material acceleration and is continuous
// across the interface.  Identical to the parent wherever the density is
// constant.  Default false = parent path.
const bool cellCorrectionAcceleration
(
    pimple.dict().getOrDefault<Switch>("cellCorrectionAcceleration", false)
);
if (cellCorrectionAcceleration)
{
    Info<< "jump: cell correction averages the face acceleration"
        << " (cellCorrectionAcceleration true)" << nl;
}

''' + anchor, "createFields.H step 3 switch")
open(p, "w").write(s)

# --------------------------------------------------------------- pEqn.H
p = os.path.join(D, "pEqn.H")
s = open(p).read()
s = sub(s,
        "    const surfaceScalarField hf(rhof*gf);\n"
        "    faceToCellAvg(hf, gC);\n"
        "    {\n"
        "        vectorField& gcf = gC.primitiveFieldRef();\n"
        "        const scalarField& rhoc = rho.primitiveField();\n"
        "        forAll(gcf, celli)\n"
        "        {\n"
        "            gcf[celli] /= rhoc[celli];\n"
        "        }\n"
        "    }",
        "    const surfaceScalarField hf(rhof*gf);\n"
        "    if (cellCorrectionAcceleration)\n"
        "    {\n"
        "        // the material acceleration is continuous across the\n"
        "        // interface, the force density is not\n"
        "        faceToCellAvg(gf, gC);\n"
        "    }\n"
        "    else\n"
        "    {\n"
        "        faceToCellAvg(hf, gC);\n"
        "        vectorField& gcf = gC.primitiveFieldRef();\n"
        "        const scalarField& rhoc = rho.primitiveField();\n"
        "        forAll(gcf, celli)\n"
        "        {\n"
        "            gcf[celli] /= rhoc[celli];\n"
        "        }\n"
        "    }",
        "pEqn.H cell correction")
open(p, "w").write(s)

print("step 3 applied from pristine, all replacements asserted")
