#!/usr/bin/env python3
"""Step 4 of interWYConsistentCenteredJumpFoam.

Diagnostic only, no change to the numerics.  Reports the residual face
acceleration

    gf = af - snGrad(p_rgh)/rho_f = (phi/|Sf| - interp(U).n)/dt ,

which is the part of the face force the exactly projected flux does not
cancel, and which the cell correction passes to the cells attenuated by
rho_f/rho_cell.  Reported as dt*gf, a velocity increment directly
comparable with maxU, split into all faces and faces touching a mixed
cell.

Runs applyJump3.py first, so all steps are regenerated from pristine.
"""

import os
import runpy

D = ("/work1/u4485952/OpenFOAM/u4485952-v2506/geometricVoFCartesian/"
     "applications/solvers/interWYConsistentCenteredJumpFoam")

runpy.run_path(os.path.join(D, "applyJump3.py"))


def sub(text, old, new, where, count=1):
    n = text.count(old)
    assert n == count, f"{n} occurrences of the anchor in {where}, expected {count}"
    return text.replace(old, new)


p = os.path.join(D, "pEqn.H")
s = open(p).read()
s = sub(s,
        "    const surfaceScalarField hf(rhof*gf);",
        "    const surfaceScalarField hf(rhof*gf);\n"
        '    #include "residualReport.H"',
        "pEqn.H residual report hook")
open(p, "w").write(s)

open(os.path.join(D, "residualReport.H"), "w").write(r'''// residual face acceleration left by the projection, as a velocity
// increment dt*gf, over all internal faces and over faces touching a
// mixed cell.  Diagnostic only.
if (jumpReport && (runTime.timeIndex() % jumpReportInterval == 0))
{
    const labelUList& own = mesh.owner();
    const labelUList& nei = mesh.neighbour();
    const scalarField& a = alpha1.primitiveField();
    const scalar aTol = 1e-6;

    scalar sAll = 0, sMix = 0, gMaxMix = 0;
    label nAll = 0, nMix = 0;

    forAll(own, facei)
    {
        const scalar dv = runTime.deltaTValue()*gf[facei];
        sAll += dv*dv;
        ++nAll;

        const bool mixed =
            (a[own[facei]] > aTol && a[own[facei]] < 1 - aTol)
         || (a[nei[facei]] > aTol && a[nei[facei]] < 1 - aTol);

        if (mixed)
        {
            sMix += dv*dv;
            gMaxMix = max(gMaxMix, mag(dv));
            ++nMix;
        }
    }

    reduce(sAll, sumOp<scalar>());
    reduce(sMix, sumOp<scalar>());
    reduce(nAll, sumOp<label>());
    reduce(nMix, sumOp<label>());
    reduce(gMaxMix, maxOp<scalar>());

    Info<< "jump: residual dt*gf  rms all "
        << (nAll ? Foam::sqrt(sAll/nAll) : 0.0)
        << "  rms mixed " << (nMix ? Foam::sqrt(sMix/nMix) : 0.0)
        << "  max mixed " << gMaxMix
        << "  on " << nMix << " faces, maxU " << gMax(mag(U)()) << endl;
}
''')

print("step 4 applied from pristine, all replacements asserted")
