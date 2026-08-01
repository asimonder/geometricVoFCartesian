#!/usr/bin/env python3
# Adds the air-side face forcing to the Dev solver. Anchored literal
# replacements; any mismatch aborts before writing.
import sys, pathlib
d = pathlib.Path('.').resolve()

def patch(fname, old, new, tag):
    p = d/fname
    src = p.read_text()
    n = src.count(old)
    if n != 1:
        sys.exit(f"PATCH ABORT: '{tag}' matched {n} times in {fname}, need 1")
    (d/(fname + '.preAirForce')).write_text(src)
    p.write_text(src.replace(old, new))
    print(f"patched {fname}: {tag}")

patch(
    'createFields.H',
    '// surface tension: active when sigma > 0; curvature via curvatureDict',
    '''// air-side driving force (wind forcing) delivered through the
// balanced face acceleration in pEqn; zero vector = inactive
const dimensionedVector airForce
(
    "airForce", dimensionSet(1, -2, -2, 0, 0), vector::zero,
    transportProperties
);
const bool airForcing = mag(airForce.value()) > VSMALL;
if (airForcing)
{
    Info<< "wy: air-side face forcing airForce = " << airForce.value()
        << " [N/m3]" << nl;
}

// surface tension: active when sigma > 0; curvature via curvatureDict''',
    'airForce read')

patch(
    'pEqn.H',
    '''    if (stActive)
    {
        af += interfaceForcePtr->surfaceTensionForce()/rhof;
    }''',
    '''    if (stActive)
    {
        af += interfaceForcePtr->surfaceTensionForce()/rhof;
    }
    if (airForcing)
    {
        // air-partial force density at faces, same weights as gravity
        // and surface tension; the cell correction gC returns it to
        // cells. Force density (1-alpha_f)*G is the air-side stress
        // budget on every face including mixed ones (the Basilisk
        // original adds a bare kinematic acceleration instead and
        // over-delivers at mixed faces by the local density ratio).
        af += (scalar(1) - fvc::interpolate(alpha1))
             *((airForce & mesh.Sf())/mesh.magSf())/rhof;
    }''',
    'af forcing')
