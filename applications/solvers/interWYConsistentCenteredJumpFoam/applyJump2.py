#!/usr/bin/env python3
"""Step 2 of interWYConsistentCenteredJumpFoam.

Face density in the pressure equation taken from the density the sweeps
actually advect (Zuzio et al. 2020, Eqs. 24-26; Kuhn et al. 2023, Eq. 30),

    alpha_f      = phi_alpha / phi        clipped to [0, 1]
    rho_f^adv    = rho1 alpha_f + rho2 (1 - alpha_f)

blended onto the linear face interpolation of alpha where phi passes
through zero, with weight w = phi^2/(phi^2 + phiRef^2).

Runs applyJump1.py first, so both steps are regenerated from pristine.
"""

import os
import runpy

D = ("/work1/u4485952/OpenFOAM/u4485952-v2506/geometricVoFCartesian/"
     "applications/solvers/interWYConsistentCenteredJumpFoam")

runpy.run_path(os.path.join(D, "applyJump1.py"))


def sub(text, old, new, where, count=1):
    n = text.count(old)
    assert n == count, f"{n} occurrences of the anchor in {where}, expected {count}"
    return text.replace(old, new)


# ------------------------------------------------------- createFields.H
p = os.path.join(D, "createFields.H")
s = open(p).read()
anchor = "autoPtr<reconstructedDistanceFunction> rdfPtr;"
s = sub(s, anchor, r'''// Face density from the density the sweeps advect.  The geometric alpha
// flux defines a face volume fraction alpha_f = phi_alpha/phi, and the
// advected face density is rho1 alpha_f + rho2 (1 - alpha_f).  This is the
// quantity the consistency argument puts in the Poisson coefficient
// (Zuzio et al. 2020, Eqs. 24-26; Kuhn et al. 2023, Eq. 30), in place of
// the parent's harmonic average of the mixture cell densities.
// Default false = parent path.
const bool advectedFaceDensity
(
    pimple.dict().getOrDefault<Switch>("advectedFaceDensity", false)
);

// blend width where phi passes through zero, as a fraction of max|phi|
const scalar advectedPhiTol
(
    pimple.dict().getOrDefault<scalar>("advectedPhiTol", 1e-3)
);

if (advectedFaceDensity && jumpFaceDensity)
{
    FatalErrorInFunction
        << "advectedFaceDensity and jumpFaceDensity are alternative face"
        << " densities; enable at most one" << exit(FatalError);
}
if (advectedFaceDensity)
{
    Info<< "jump: face density from the advected alpha flux"
        << " (advectedFaceDensity true, advectedPhiTol "
        << advectedPhiTol << ")" << nl;
}

''' + anchor, "createFields.H step 2 switches")
open(p, "w").write(s)

# --------------------------------------------------------------- pEqn.H
p = os.path.join(D, "pEqn.H")
s = open(p).read()
s = sub(s,
        "    if (jumpFaceDensity)\n"
        "    {\n"
        '        #include "jumpFaceDensity.H"\n'
        "    }",
        "    if (jumpFaceDensity)\n"
        "    {\n"
        '        #include "jumpFaceDensity.H"\n'
        "    }\n"
        "    if (advectedFaceDensity)\n"
        "    {\n"
        '        #include "advectedFaceDensity.H"\n'
        "    }",
        "pEqn.H step 2 hook")
open(p, "w").write(s)

# ------------------------------------------------- advectedFaceDensity.H
open(os.path.join(D, "advectedFaceDensity.H"), "w").write(r'''// Face density from the density the sweeps advect.
//
//     alpha_f   = phi_alpha/phi, clipped to [0, 1]
//     rho_f     = rho1 alpha_f + rho2 (1 - alpha_f)
//
// Where phi passes through zero the ratio is meaningless, so the value is
// blended onto the linear face interpolation of alpha with the weight
// w = phi^2/(phi^2 + phiRef^2), phiRef a fixed fraction of max|phi|.
{
    const surfaceScalarField& aPhi = advector.alphaPhi();

    const dimensionedScalar phiSmall("phiSmall", phi.dimensions(), VSMALL);
    const dimensionedScalar phiRef
    (
        "phiRef",
        phi.dimensions(),
        advectedPhiTol*(gMax(mag(phi.primitiveField())) + VSMALL)
    );

    surfaceScalarField alphaf
    (
        min(max(aPhi/stabilise(phi, phiSmall), scalar(0)), scalar(1))
    );

    const surfaceScalarField w(sqr(phi)/(sqr(phi) + sqr(phiRef)));
    alphaf = w*alphaf + (scalar(1) - w)*fvc::interpolate(alpha1);

    rhof = rho1*alphaf + rho2*(scalar(1) - alphaf);

    if (jumpReport && (runTime.timeIndex() % jumpReportInterval == 0))
    {
        Info<< "jump: advected face density, rho_f in ["
            << gMin(rhof.primitiveField()) << ", "
            << gMax(rhof.primitiveField()) << "] kg/m3, phiRef "
            << phiRef.value() << endl;
    }
}
''')

print("step 2 applied from pristine, all replacements asserted")
