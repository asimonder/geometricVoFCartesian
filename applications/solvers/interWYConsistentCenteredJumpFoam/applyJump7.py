#!/usr/bin/env python3
"""Step 7 of interWYConsistentCenteredJumpFoam.

Queutey and Visonneau (2007) hydrostatic-type discontinuity treatment.

Their jump conditions at the interface Gamma are Eqs. (34a) and (34b),

    [p] = 0,        [grad p / rho] = 0,

and the normal gradient over density at the face is their Eqs. (37), (38),

    (grad p . n / rho)_f = (1/rho_hat) (p_R - p_L)/h ,
    rho_hat = (h^- rho^- + h^+ rho^+)/h                          Eq. (36b)

which on a uniform orthogonal mesh is the arithmetic mean of the two side
densities.  Their Sec. 4.2.1 places the discontinuity on a mesh face with
pure phases on either side, so at that face rho_hat = (rho1 + rho2)/2.

Translated into the p_rgh variable used here, with p = p_rgh + rho gh and
snGrad(gh) = g.n on a uniform mesh,

    snGrad(rho gh) = gh_f snGrad(rho) + rho_arith,f (g.n)

so their face acceleration g.n - snGrad(p)/rho_hat reduces to
-gh_f snGrad(rho)/rho_hat, the g.n cancelling identically.  Their Eq. (36a)
face pressure value has no counterpart in a projection cycle, where the
pressure enters only through the face normal gradient.

The implementation therefore selects the face density at the faces the
interface crosses, reusing the distance-function machinery already tested
in step 1, through a new word

    jumpFaceDensityModel  oneSided | queutey

with oneSided the default, so the existing switch is unchanged.

Runs applyJump6.py first, so all steps are regenerated from pristine.
"""

import os
import runpy

D = ("/work1/u4485952/OpenFOAM/u4485952-v2506/geometricVoFCartesian/"
     "applications/solvers/interWYConsistentCenteredJumpFoam")

runpy.run_path(os.path.join(D, "applyJump6.py"))


def sub(text, old, new, where, count=1):
    n = text.count(old)
    assert n == count, f"{n} occurrences of the anchor in {where}, expected {count}"
    return text.replace(old, new)


# ------------------------------------------------------- createFields.H
p = os.path.join(D, "createFields.H")
s = open(p).read()
anchor = "autoPtr<reconstructedDistanceFunction> rdfPtr;"
s = sub(s, anchor, r'''// which face density is used at the faces the interface crosses.
//   oneSided : lambda/rho_s(P) + (1-lambda)/rho_s(N), the 1D-exact series
//              value (Vukcevic et al. 2017, Eqs. 39-49)
//   queutey  : the arithmetic mean of the two side densities, rho_hat of
//              Queutey and Visonneau (2007), Eqs. (36b), (37), (38)
const word jumpFaceDensityModel
(
    pimple.dict().getOrDefault<word>("jumpFaceDensityModel", "oneSided")
);
if (jumpFaceDensity)
{
    if (jumpFaceDensityModel != "oneSided"
     && jumpFaceDensityModel != "queutey")
    {
        FatalErrorInFunction
            << "jumpFaceDensityModel must be oneSided or queutey, got "
            << jumpFaceDensityModel << exit(FatalError);
    }
    Info<< "jump: interface face density model " << jumpFaceDensityModel
        << nl;
}

''' + anchor, "createFields.H step 7 word")
open(p, "w").write(s)

# ----------------------------------------------------- jumpFaceDensity.H
p = os.path.join(D, "jumpFaceDensity.H")
s = open(p).read()
s = sub(s,
        "    auto oneSided = [r1, r2](const scalar psiP, const scalar psiN)\n"
        "    {\n"
        "        const scalar den = psiP - psiN;\n"
        "        const scalar lambda = min(max(psiP/den, scalar(0)), scalar(1));\n"
        "        const scalar rP = (psiP > 0) ? r1 : r2;\n"
        "        const scalar rN = (psiN > 0) ? r1 : r2;\n"
        "        return scalar(1)/(lambda/rP + (scalar(1) - lambda)/rN);\n"
        "    };",
        "    const bool qv = (jumpFaceDensityModel == \"queutey\");\n"
        "    auto oneSided = [r1, r2, qv](const scalar psiP, const scalar psiN)\n"
        "    {\n"
        "        const scalar rP = (psiP > 0) ? r1 : r2;\n"
        "        const scalar rN = (psiN > 0) ? r1 : r2;\n"
        "        if (qv)\n"
        "        {\n"
        "            // Queutey and Visonneau (2007) Eq. (36b) on a uniform\n"
        "            // orthogonal mesh with pure phases either side\n"
        "            return 0.5*(rP + rN);\n"
        "        }\n"
        "        const scalar den = psiP - psiN;\n"
        "        const scalar lambda = min(max(psiP/den, scalar(0)), scalar(1));\n"
        "        return scalar(1)/(lambda/rP + (scalar(1) - lambda)/rN);\n"
        "    };",
        "jumpFaceDensity.H model selection")
open(p, "w").write(s)

print("step 7 applied from pristine, all replacements asserted")
