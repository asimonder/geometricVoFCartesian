/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2026 Asim Onder
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    interIsoConsistentCenteredFoam

Description
    Mass-momentum consistent two-phase solver on collocated meshes,
    general-mesh companion of interWYConsistentCenteredFoam:

      alphaEqn:  volume fraction advected by isoAdvector (Roenby,
                 Bredmose, Jasak 2016), staggered dt/2 once at the dt
                 plateau; the mixture momentum is advected explicitly
                 with the SAME geometric phase volume fluxes (MCVOF,
                 Zhang, Popinet, Ling 2020, Eqs. 12-19; unsplit
                 consistency per Kuhn, Deskos, Sprague 2023), donor
                 values from per-phase masked slopes (centered
                 unlimited by default) with BCG-lite half-step tracing.
      UEqn:      u-star = q-star over rho(alpha_new), recovered in the
                 advection; implicit viscous step on u-star + dt g
                 (correction form).
      pEqn:      face flux from the cell velocity plus the balanced
                 face forces (gravity, surface tension, air forcing);
                 ONE exact face projection with HARMONIC face density
                 rAUf = dt/rho_f, rho_f = 1/interp(1/rho); no
                 Rhie-Chow, no ddtCorr, no correctors. Cell velocities
                 corrected by the approximate projection: the
                 weight-free face force density gathered to cells by
                 fvc::reconstruct (general meshes; reduces to the
                 per-direction two-face average on uniform Cartesian
                 grids) or by the explicit axis average
                 (cellCorrection average), divided by the CELL density.

    Surface tension uses interfaceForces (curvatureDict), active when
    sigma > 0; sigma = 0 skips construction entirely.

    nOuterCorrectors and nAlphaSubCycles other than 1 are rejected:
    the consistency of the cycle is structural.

    Topology change: the correctPhi.H flux restoration of the centered
    family (single global reference cell) is retained; the MCVOF
    momentum advection is untested under mesh change and warns.
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "isoAdvection.H"
#include "immiscibleIncompressibleTwoPhaseGenericMixture.H"
#include "pimpleControl.H"
#include "interfaceForces.H"
#include "CorrectPhi.H"
#include "adjustPhi.H"
#include "fixedValueFvPatchFields.H"
#include "mathematicalConstants.H"
#include "emptyFvPatch.H"
#include "fvOptions.H"
#include "incompressibleInterPhaseTransportModel.H"
#include "mcvofSlopes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Mass-momentum consistent two-phase solver on collocated meshes"
        " (isoAdvector + MCVOF consistent transport, centered"
        " single-projection cycle, harmonic face-density weights)."
    );

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "initContinuityErrs.H"
    #include "createDyMControls.H"
    #include "createFields.H"

    #include "CourantNo.H"
    #include "setInitialDeltaT.H"

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readDyMControls.H"
        if (stActive) { maxDeltaT = min(maxDeltaT, dtCap); }
        #include "CourantNo.H"
        #include "alphaCourantNo.H"
        #include "setDeltaT.H"

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        const dimensionedScalar dtDim("dt", dimTime, runTime.deltaTValue());

        // one-time dt/2 alpha staggering at the dt plateau (the
        // staggered-in-time volume fraction of the centered cycle)
        if (staggerAlpha && !staggerDone)
        {
            const scalar dtNow = runTime.deltaTValue();
            if (dtPrev > 0 && mag(dtNow - dtPrev) < 1e-3*dtNow)
            {
                Info<< "isoCentered: dt plateau; staggering alpha by dt/2 = "
                    << 0.5*dtNow << " s" << nl;
                phi *= 0.5;
                {
                    #include "alphaSuSp.H"
                    advector.advect(Sp, Su);
                }
                phi *= 2.0;
                alpha2 = 1.0 - alpha1;
                alpha2.correctBoundaryConditions();
                mixture.correct();
                rho == alpha1*rho1 + alpha2*rho2;
                staggerDone = true;
            }
            dtPrev = dtNow;
        }

        // --- Pressure-velocity corrector loop (single pass enforced)
        while (pimple.loop())
        {
            if (pimple.firstIter() || moveMeshOuterCorrectors)
            {
                if (mesh.dynamic())
                {
                    // store the PLIC planes so a topology change can
                    // rebuild alpha sharply
                    advector.surf().reconstruct();
                }

                mesh.update();

                if (mesh.changing())
                {
                    gh = (g & mesh.C()) - ghRef;
                    ghf = (g & mesh.Cf()) - ghRef;

                    advector.surf().mapAlphaField();
                    alpha2 = 1.0 - alpha1;
                    alpha2.correctBoundaryConditions();
                    rho == alpha1*rho1 + alpha2*rho2;
                    rho.correctBoundaryConditions();
                    rho.oldTime() = rho;
                    alpha2.oldTime() = alpha2;

                    if (correctPhi)
                    {
                        // single-global-reference projection; see
                        // correctPhi.H
                        #include "correctPhi.H"
                        #include "continuityErrs.H"
                    }
                }
            }

            #include "alphaControls.H"
            #include "alphaEqnSubCycle.H"

            mixture.correct();

            if (stActive)
            {
                interfaceForcePtr->correct();
            }

            #include "UEqn.H"
            #include "pEqn.H"

            if (turbulenceActive && pimple.turbCorr())
            {
                turbulence->correct();
            }
        }

        runTime.write();

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;
    return 0;
}


// ************************************************************************* //
