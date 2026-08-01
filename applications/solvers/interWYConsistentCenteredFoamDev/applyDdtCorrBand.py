#!/usr/bin/env python3
# Extends the ddtCorr variant with a scale coefficient and an
# interface-band restriction: wyFluxDdtCorrCoeff (default 1.0),
# wyFluxDdtCorrBand (default false; term applied only on faces
# adjacent to mixed cells). Run after applyDdtCorr.py.
import sys, pathlib
d = pathlib.Path('.').resolve()

def patch(fname, old, new, tag):
    p = d/fname
    src = p.read_text()
    n = src.count(old)
    if n != 1:
        sys.exit(f"PATCH ABORT: '{tag}' matched {n} times in {fname}, need 1")
    (d/(fname + '.preBand2')).write_text(src)
    p.write_text(src.replace(old, new))
    print(f"patched {fname}: {tag}")

patch(
    'createFields.H',
    'if (fluxDdtCorr)\n{\n    Info<< "wy: phiHbyA carries the ddtCorr face-flux consistency"\n        << " term (variant, wyFluxDdtCorr true)" << nl;\n}',
    '''const scalar ddtCorrCoeff
(
    pimple.dict().getOrDefault<scalar>("wyFluxDdtCorrCoeff", 1.0)
);
const bool fluxDdtCorrBand
(
    pimple.dict().getOrDefault<Switch>("wyFluxDdtCorrBand", false)
);
if (fluxDdtCorr)
{
    Info<< "wy: phiHbyA carries the ddtCorr face-flux consistency"
        << " term (variant, wyFluxDdtCorr true, coeff "
        << ddtCorrCoeff
        << (fluxDdtCorrBand ? ", interface band only)" : ", all faces)")
        << nl;
}''',
    'coeff and band flags')

patch(
    'pEqn.H',
    '''    if (fluxDdtCorr)
    {
        phiHbyA += dtDim*fvc::ddtCorr(U, phi);
    }''',
    '''    if (fluxDdtCorr)
    {
        const surfaceScalarField corrF(dtDim*fvc::ddtCorr(U, phi));
        if (fluxDdtCorrBand)
        {
            // flux memory only on faces adjacent to mixed cells,
            // where the misdeposit is generated
            const volScalarField bandCell
            (
                "wyBandCell",
                pos(alpha1*(scalar(1) - alpha1) - SMALL)
            );
            phiHbyA +=
                ddtCorrCoeff
               *pos(fvc::interpolate(bandCell) - SMALL)*corrF;
        }
        else
        {
            phiHbyA += ddtCorrCoeff*corrF;
        }
    }''',
    'banded term')
