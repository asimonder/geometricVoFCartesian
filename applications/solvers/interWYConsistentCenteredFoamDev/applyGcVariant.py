#!/usr/bin/env python3
# Adds a runtime-selectable cell-correction variant to the Dev solver:
# wyCellCorrectionReconstruct true replaces the axis-averaged face-force
# gather with fvc::reconstruct (the old conserving solver's operator).
# Default false keeps the production path bit for bit. Run after the
# other apply*.py scripts. Anchored replacements; abort on mismatch.
import sys, pathlib
d = pathlib.Path('.').resolve()

def patch(fname, old, new, tag):
    p = d/fname
    src = p.read_text()
    n = src.count(old)
    if n != 1:
        sys.exit(f"PATCH ABORT: '{tag}' matched {n} times in {fname}, need 1")
    (d/(fname + '.preGcVar')).write_text(src)
    p.write_text(src.replace(old, new))
    print(f"patched {fname}: {tag}")

patch(
    'createFields.H',
    '// turbulence transport (interIsoCartFoam pattern): constructed only',
    '''// cell-correction variant: reconstruct-based face-force gather
// (the old conserving solver's operator) instead of the axis average;
// default false = production path bit for bit
const bool gcReconstruct
(
    pimple.dict().getOrDefault<Switch>("wyCellCorrectionReconstruct", false)
);
if (gcReconstruct)
{
    Info<< "wy: cell correction gC via fvc::reconstruct"
        << " (variant, wyCellCorrectionReconstruct true)" << nl;
}

// turbulence transport (interIsoCartFoam pattern): constructed only''',
    'gcReconstruct flag')

patch(
    'pEqn.H',
    '''    faceToCellAvg(hf, gC);
    {''',
    '''    if (gcReconstruct)
    {
        gC.primitiveFieldRef() =
            fvc::reconstruct(hf*mesh.magSf())().primitiveField();
    }
    else
    {
        faceToCellAvg(hf, gC);
    }
    {''',
    'gC variant branch')
