#!/usr/bin/env python3
"""Step 8 of interWYConsistentCenteredJumpFoam.

Second-order one-sided cell correction.

Step 6 excluded the face touching a mixed cell and gave the remaining face
the whole weight.  That value sits half a cell from the cell centre, so it
is first order where the parent's average is second order, in every cell
beside the band on both sides.  This step uses the one-sided second-order
difference from the two clean faces on the same side.  With v1 the adjacent
clean face value and v2 the next one out, extrapolation to the cell centre
gives

    g_c = (3/2) v1 - (1/2) v2 ,

and since the neighbour across the clean face stores
S = (v1 + v2)/2 when both of its faces in that direction are clean,

    g_c = 2 v1 - S .

Checked against a linear face field: v(x) = x, faces at h/2 and 3h/2 give
g_c = 0 at the cell centre, from either side.

Where the second face is missing, or the neighbour does not have both of
its faces clean, the first-order value is kept.  Where both faces of a
direction are clean, or both touch a mixed cell, the parent average stands.

Selected by

    cellCorrectionOneSidedOrder   1 | 2      (default 1, the step-6 form)

Runs applyJump7.py first, so all steps are regenerated from pristine.
"""

import os
import runpy

D = ("/work1/u4485952/OpenFOAM/u4485952-v2506/geometricVoFCartesian/"
     "applications/solvers/interWYConsistentCenteredJumpFoam")

runpy.run_path(os.path.join(D, "applyJump7.py"))


def sub(text, old, new, where, count=1):
    n = text.count(old)
    assert n == count, f"{n} occurrences of the anchor in {where}, expected {count}"
    return text.replace(old, new)


p = os.path.join(D, "createFields.H")
s = open(p).read()

s = sub(s,
        "const scalar cellCorrectionAlphaTol\n"
        "(\n"
        '    pimple.dict().getOrDefault<scalar>("cellCorrectionAlphaTol", 1e-6)\n'
        ");",
        "const scalar cellCorrectionAlphaTol\n"
        "(\n"
        '    pimple.dict().getOrDefault<scalar>("cellCorrectionAlphaTol", 1e-6)\n'
        ");\n"
        "// 1 = adjacent clean face at full weight (first order)\n"
        "// 2 = one-sided second-order difference from the two clean faces\n"
        "const label cellCorrectionOneSidedOrder\n"
        "(\n"
        '    pimple.dict().getOrDefault<label>("cellCorrectionOneSidedOrder", 1)\n'
        ");",
        "createFields.H order control")

s = sub(s,
        "[&mesh, &alpha1, &cellCorrectionAlphaTol]",
        "[&mesh, &alpha1, &cellCorrectionAlphaTol, &cellCorrectionOneSidedOrder]",
        "createFields.H lambda capture")

s = sub(s,
        """    vectorField cleanSum(gc.size(), Zero);
    vectorField nAll(gc.size(), Zero);
    vectorField nClean(gc.size(), Zero);""",
        """    vectorField cleanSum(gc.size(), Zero);
    vectorField nAll(gc.size(), Zero);
    vectorField nClean(gc.size(), Zero);
    // cell across the clean face of each direction, for the second-order
    // one-sided form; only read where a direction has exactly one clean
    // face, so a cell with two clean faces overwriting it is harmless
    labelList cleanNbr(3*gc.size(), -1);""",
        "createFields.H one-sided storage")

s = sub(s,
        """            const scalar v = gfld[facei]*sign(nHat[d]);
            add(own[facei], d, v, clean);
            add(nei[facei], d, v, clean);""",
        """            const scalar v = gfld[facei]*sign(nHat[d]);
            add(own[facei], d, v, clean);
            add(nei[facei], d, v, clean);
            if (clean)
            {
                cleanNbr[3*own[facei] + d] = nei[facei];
                cleanNbr[3*nei[facei] + d] = own[facei];
            }""",
        "createFields.H clean neighbour")

s = sub(s,
        """    forAll(gc, celli)
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
    gcell.correctBoundaryConditions();""",
        """    forAll(gc, celli)
    {
        for (direction d = 0; d < 3; ++d)
        {
            const scalar na = nAll[celli][d];
            const scalar nc = nClean[celli][d];
            if (!(na > 0 && nc > 0 && nc < na)) continue;

            // v1, the adjacent clean face value at full weight
            const scalar v1 = cleanSum[celli][d]*(na/nc);
            gc[celli][d] = v1;

            if (cellCorrectionOneSidedOrder < 2) continue;

            const label nb = cleanNbr[3*celli + d];
            if (nb < 0 || nClean[nb][d] < 1.5) continue;

            // the neighbour holds S = (v1 + v2)/2 with both faces clean
            gc[celli][d] = 2.0*v1 - cleanSum[nb][d];
        }
    }
    gcell.correctBoundaryConditions();""",
        "createFields.H one-sided finalisation")

open(p, "w").write(s)
print("step 8 applied from pristine, all replacements asserted")
