#!/usr/bin/env python3
# Adds wyFluxDdtCorr (default false = production path): the stock
# ddtCorr face-flux consistency term in phiHbyA, restoring the memory
# of the previous projected flux that the per-step re-interpolation
# of the cell velocity discards. Run after the other apply*.py.
import sys, pathlib
d = pathlib.Path('.').resolve()

def patch(fname, old, new, tag):
    p = d/fname
    src = p.read_text()
    n = src.count(old)
    if n != 1:
        sys.exit(f"PATCH ABORT: '{tag}' matched {n} times in {fname}, need 1")
    (d/(fname + '.preDdtCorr')).write_text(src)
    p.write_text(src.replace(old, new))
    print(f"patched {fname}: {tag}")

patch(
    'createFields.H',
    '// pressure-pass count: additional momentum-pressure coupling',
    '''// face-flux consistency term (stock ddtCorr) in phiHbyA;
// default false = production path
const bool fluxDdtCorr
(
    pimple.dict().getOrDefault<Switch>("wyFluxDdtCorr", false)
);
if (fluxDdtCorr)
{
    Info<< "wy: phiHbyA carries the ddtCorr face-flux consistency"
        << " term (variant, wyFluxDdtCorr true)" << nl;
}

// pressure-pass count: additional momentum-pressure coupling''',
    'ddtCorr flag')

patch(
    'pEqn.H',
    '''    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        (fvc::interpolate(U) & mesh.Sf()) + dtDim*af*mesh.magSf()
    );''',
    '''    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        (fvc::interpolate(U) & mesh.Sf()) + dtDim*af*mesh.magSf()
    );
    if (fluxDdtCorr)
    {
        phiHbyA += dtDim*fvc::ddtCorr(U, phi);
    }''',
    'ddtCorr term')
