#!/usr/bin/env python3
# Adds turbulence transport (interIsoCartFoam pattern) with a laminar
# bypass: cases selecting simulationType laminar keep the original
# viscous path bit for bit; a real model routes through the stock
# incompressibleInterPhaseTransportModel. Run AFTER applyDevPatch.py
# and applyAirForce.py. Anchored literal replacements; abort on
# mismatch.
import sys, pathlib
d = pathlib.Path('.').resolve()

def patch(fname, old, new, tag, backup=True):
    p = d/fname
    src = p.read_text()
    n = src.count(old)
    if n != 1:
        sys.exit(f"PATCH ABORT: '{tag}' matched {n} times in {fname}, need 1")
    if backup:
        (d/(fname + '.preTurb')).write_text(src)
    p.write_text(src.replace(old, new))
    print(f"patched {fname}: {tag}")

# ---- createFields.H: guarded construction ------------------------------
patch(
    'createFields.H',
    '''if (airForcing)
{
    Info<< "wy: air-side face forcing airForce = " << airForce.value()
        << " [N/m3]" << nl;
}''',
    '''if (airForcing)
{
    Info<< "wy: air-side face forcing airForce = " << airForce.value()
        << " [N/m3]" << nl;
}

// turbulence transport (interIsoCartFoam pattern): constructed only
// when turbulenceProperties selects a model; laminar cases keep the
// original viscous path bit for bit (laminar bypass)
bool turbulenceActive = false;
{
    IOobject turbIO
    (
        "turbulenceProperties", runTime.constant(), mesh,
        IOobject::MUST_READ, IOobject::NO_WRITE, false
    );
    if (turbIO.typeHeaderOk<IOdictionary>(true))
    {
        IOdictionary turbProps(turbIO);
        turbulenceActive =
            turbProps.getOrDefault<word>("simulationType", "laminar")
         != "laminar";
    }
}

typedef incompressibleInterPhaseTransportModel
    <
        immiscibleIncompressibleTwoPhaseGenericMixture
    > transportModelType;

autoPtr<surfaceScalarField> rhoPhiPtr;
autoPtr<transportModelType> turbulence;
if (turbulenceActive)
{
    rhoPhiPtr.reset
    (
        new surfaceScalarField
        (
            IOobject
            (
                "rhoPhi", runTime.timeName(), mesh,
                IOobject::NO_READ, IOobject::NO_WRITE
            ),
            fvc::interpolate(rho)*phi
        )
    );
    turbulence.reset
    (
        new transportModelType(rho, U, phi, rhoPhiPtr(), mixture)
    );
    Info<< "wy: turbulence transport ACTIVE"
        << " (mass-consistent rhoPhi from the geometric alpha flux)"
        << nl;
}''',
    'turbulence construction')

# ---- UEqn.H: laminar bypass branch --------------------------------------
patch(
    'UEqn.H',
    '''    U += dtDim*gC;
    {
        const surfaceScalarField muf(mixture.muf());''',
    '''    U += dtDim*gC;
    if (turbulenceActive)
    {
        // stock stress operator; the laminar branch below is the
        // original path and stays digit-identical when no model is
        // selected
        volVectorField Uprov("Uprov", U);
        fvVectorMatrix UEqn
        (
            fvm::Sp(rho/dtDim, U)
          + turbulence->divDevRhoReff(rho, U)
         ==
            rho/dtDim*Uprov
        );
        UEqn -= fvOptions(rho, U);
        fvOptions.constrain(UEqn);
        UEqn.solve();
    }
    else
    {
        const surfaceScalarField muf(mixture.muf());''',
    'UEqn bypass')

# ---- alphaEqn.H: mass-consistent flux update ----------------------------
patch(
    'alphaEqn.H',
    '''    alpha2 = 1.0 - alpha1;
    rho == alpha1*rho1 + alpha2*rho2;

    Info<< "Phase-1 volume fraction = "''',
    '''    alpha2 = 1.0 - alpha1;
    rho == alpha1*rho1 + alpha2*rho2;

    if (turbulenceActive)
    {
        // mass-consistent face flux from the geometric alpha flux
        const surfaceScalarField& aPhi = advector.alphaPhi();
        rhoPhiPtr() = rho1*aPhi + rho2*(phi - aPhi);
    }

    Info<< "Phase-1 volume fraction = "''',
    'rhoPhi update')

# ---- main .C (already generated): include + correct() -------------------
patch(
    'interWYConsistentCenteredFoamDev.C',
    '#include "fvOptions.H"',
    '#include "fvOptions.H"\n#include "incompressibleInterPhaseTransportModel.H"',
    'transport include', backup=False)

patch(
    'interWYConsistentCenteredFoamDev.C',
    '''            #include "UEqn.H"
            #include "pEqn.H"''',
    '''            #include "UEqn.H"
            #include "pEqn.H"

            if (turbulenceActive && pimple.turbCorr())
            {
                turbulence->correct();
            }''',
    'turbulence correct', backup=False)
