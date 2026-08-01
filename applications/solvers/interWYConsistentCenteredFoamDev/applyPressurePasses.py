#!/usr/bin/env python3
# Adds wyPressurePasses (default 1 = production path): additional
# momentum-pressure coupling passes inside the centered cycle. Each
# pass rebuilds phiHbyA from the corrected velocity and re-projects;
# the cell correction REPLACES the previous pass's (U = UPreCorr +
# dt*gC), so passes converge the coupling instead of double-stepping.
# Run after the other apply*.py scripts. Anchored; abort on mismatch.
import sys, pathlib
d = pathlib.Path('.').resolve()

def patch(fname, old, new, tag):
    p = d/fname
    src = p.read_text()
    n = src.count(old)
    if n != 1:
        sys.exit(f"PATCH ABORT: '{tag}' matched {n} times in {fname}, need 1")
    (d/(fname + '.prePasses')).write_text(src)
    p.write_text(src.replace(old, new))
    print(f"patched {fname}: {tag}")

patch(
    'createFields.H',
    '// cell-correction variant: reconstruct-based face-force gather',
    '''// pressure-pass count: additional momentum-pressure coupling
// passes per step (PISO-corrector analogue); 1 = production path
const label wyPressurePasses
(
    pimple.dict().getOrDefault<label>("wyPressurePasses", 1)
);
if (wyPressurePasses < 1)
{
    FatalErrorInFunction
        << "wyPressurePasses must be >= 1" << exit(FatalError);
}
if (wyPressurePasses > 1)
{
    Info<< "wy: " << wyPressurePasses
        << " momentum-pressure coupling passes per step (variant)" << nl;
}

// cell-correction variant: reconstruct-based face-force gather''',
    'passes flag')

patch(
    'pEqn.H',
    '    U += dtDim*gC;',
    '    U = UPreCorr + dtDim*gC;',
    'replacing correction')

patch(
    'interWYConsistentCenteredFoamDev.C',
    '''            #include "UEqn.H"
            #include "pEqn.H"''',
    '''            #include "UEqn.H"
            {
                const volVectorField UPreCorr("UPreCorr", U);
                for (label wyPass = 0; wyPass < wyPressurePasses; ++wyPass)
                {
                    #include "pEqn.H"
                }
            }''',
    'pass loop')
