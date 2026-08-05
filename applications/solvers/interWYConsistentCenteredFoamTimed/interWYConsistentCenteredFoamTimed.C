/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | www.openfoam.com
    \\  /    A nd           |
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
    interWYConsistentCenteredFoamTimed

Description
    Timed clone of interWYConsistentCenteredFoam. The algorithm is
    unchanged. The only additions are wall-clock timers around the
    stage includes of the time loop, and the report they print.

    The clock is a single clockTime instance. Consecutive calls to
    timeIncrement() partition the step exactly, so the stage totals sum
    to the loop total up to the cost of the calls themselves.

    Stages, in the order they are timed:
      ctl    controls, Courant numbers, setDeltaT
      mesh   mesh.update, alpha remap, correctPhi
      adv    alphaEqnSubCycle.H          (Weymouth-Yue sweeps, q1, q2)
      surf   mixture.correct, curvature and surface tension
      visc   UEqn.H                      (velocity recovery, viscous solve)
      pres   pEqn.H                      (face flux, Poisson solve, gC)
      turb   turbulence->correct
      write  runTime.write

    Per step the master rank prints one wyTimers line. After the loop the
    cumulative totals are reduced with maxOp over ranks and printed as
    seconds and as a percentage of the accounted total.

    NOTE: pres covers the whole of pEqn.H, which is the face-force
    assembly, the Poisson solve and the averaged cell correction. It is
    an upper bound on the linear-solve cost, not the linear-solve cost.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "geometricAdvector.H"
#include "splitAdvector.H"
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
#include "clockTime.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Mass-momentum consistent two-phase solver on collocated grids"
        " (Weymouth-Yue consistent transport, single-projection cycle,"
        " harmonic face-density weights). Timed build."
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

    Info<< "interWYConsistentCenteredFoam: correctPhi rebuilt with"
        << " the single global reference cell (pRefCell), conditional flux rebuild (v3); library"
        << " CorrectPhi rank-local setReference(0,0) removed" << endl;

    // --- stage timers (timed build only) -----------------------------------
    clockTime wyClk;
    clockTime wyLoopClk;
    scalar wyCtl = 0, wyMesh = 0, wyAdv = 0, wySurf = 0,
           wyVisc = 0, wyPres = 0, wyTurb = 0, wyWrite = 0;
    label wySteps = 0;
    const label wyReport =
        runTime.controlDict().getOrDefault<label>("wyTimerInterval", 1);
    Info<< "wyTimers: stage timing active, report interval "
        << wyReport << " step(s)" << endl;
    // -----------------------------------------------------------------------

    Info<< "\nStarting time loop\n" << endl;

    wyLoopClk.timeIncrement();

    while (runTime.run())
    {
        wyClk.timeIncrement();

        #include "readDyMControls.H"
        if (stActive) { maxDeltaT = min(maxDeltaT, dtCap); }
        #include "CourantNo.H"
        #include "alphaCourantNo.H"
        #include "setDeltaT.H"

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        const dimensionedScalar dtDim("dt", dimTime, runTime.deltaTValue());

        scalar sCtl = wyClk.timeIncrement();
        scalar sMesh = 0, sAdv = 0, sSurf = 0, sVisc = 0, sPres = 0, sTurb = 0;

        // --- Pressure-velocity corrector loop (single pass enforced)
        while (pimple.loop())
        {
            if (pimple.firstIter() || moveMeshOuterCorrectors)
            {
                if (mesh.dynamic())
                {
                    // store the PLIC planes so a topology change can
                    // rebuild alpha sharply (interIsoDynCartFoam pattern)
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
                        // the remapped flux is not solenoidal; the WY
                        // compression identity needs div(phi) = 0 before
                        // the advection consumes it
                        // single-global-reference projection; see correctPhi.H
                        #include "correctPhi.H"
                        #include "continuityErrs.H"
                    }
                }
            }

            sMesh += wyClk.timeIncrement();

            #include "alphaControls.H"
            #include "alphaEqnSubCycle.H"

            sAdv += wyClk.timeIncrement();

            mixture.correct();

            if (stActive)
            {
                interfaceForcePtr->correct();
            }

            sSurf += wyClk.timeIncrement();

            #include "UEqn.H"

            sVisc += wyClk.timeIncrement();

            #include "pEqn.H"

            sPres += wyClk.timeIncrement();

            if (turbulenceActive && pimple.turbCorr())
            {
                turbulence->correct();
            }

            sTurb += wyClk.timeIncrement();
        }

        runTime.write();

        const scalar sWrite = wyClk.timeIncrement();

        // --- accumulate and report ----------------------------------------
        wyCtl += sCtl;   wyMesh += sMesh;  wyAdv += sAdv;    wySurf += sSurf;
        wyVisc += sVisc; wyPres += sPres;  wyTurb += sTurb;  wyWrite += sWrite;
        ++wySteps;

        if ((wySteps % wyReport) == 0)
        {
            const scalar sSum =
                sCtl + sMesh + sAdv + sSurf + sVisc + sPres + sTurb + sWrite;
            const scalar cSum =
                wyCtl + wyMesh + wyAdv + wySurf
              + wyVisc + wyPres + wyTurb + wyWrite;

            Info<< "wyTimers  ctl " << sCtl
                << "  mesh " << sMesh
                << "  adv " << sAdv
                << "  surf " << sSurf
                << "  visc " << sVisc
                << "  pres " << sPres
                << "  turb " << sTurb
                << "  write " << sWrite
                << "  step " << sSum << " s" << endl;

            if (cSum > SMALL)
            {
                Info<< "wyTimers% adv " << 100.0*wyAdv/cSum
                    << "  visc " << 100.0*wyVisc/cSum
                    << "  pres " << 100.0*wyPres/cSum
                    << "  surf " << 100.0*wySurf/cSum
                    << "  mesh " << 100.0*wyMesh/cSum
                    << "  ctl " << 100.0*wyCtl/cSum
                    << "  turb " << 100.0*wyTurb/cSum
                    << "  write " << 100.0*wyWrite/cSum
                    << "   (cumulative, " << wySteps << " steps)" << endl;
            }
        }

        runTime.printExecutionTime(Info);
    }

    // --- final report, worst rank ------------------------------------------
    {
        const scalar wyLoop = returnReduce(wyLoopClk.elapsedTime(), maxOp<scalar>());
        const scalar mCtl   = returnReduce(wyCtl,   maxOp<scalar>());
        const scalar mMesh  = returnReduce(wyMesh,  maxOp<scalar>());
        const scalar mAdv   = returnReduce(wyAdv,   maxOp<scalar>());
        const scalar mSurf  = returnReduce(wySurf,  maxOp<scalar>());
        const scalar mVisc  = returnReduce(wyVisc,  maxOp<scalar>());
        const scalar mPres  = returnReduce(wyPres,  maxOp<scalar>());
        const scalar mTurb  = returnReduce(wyTurb,  maxOp<scalar>());
        const scalar mWrite = returnReduce(wyWrite, maxOp<scalar>());
        const scalar acc =
            mCtl + mMesh + mAdv + mSurf + mVisc + mPres + mTurb + mWrite;
        const scalar den = max(acc, SMALL);
        const scalar nrm = max(scalar(wySteps), scalar(1));

        Info<< nl << "wyTimers final report (max over ranks)" << nl
            << "    steps           " << wySteps << nl
            << "    stage            seconds     percent    s/step" << nl;

        #define WY_ROW(name, val)                                             \
            Info<< "    " << name << "   " << (val)                           \
                << "   " << 100.0*(val)/den                                   \
                << "   " << (val)/nrm << nl;

        WY_ROW("advection      ", mAdv)
        WY_ROW("pressure       ", mPres)
        WY_ROW("viscous        ", mVisc)
        WY_ROW("surfaceTension ", mSurf)
        WY_ROW("meshUpdate     ", mMesh)
        WY_ROW("control        ", mCtl)
        WY_ROW("turbulence     ", mTurb)
        WY_ROW("write          ", mWrite)
        #undef WY_ROW

        Info<< "    accounted       " << acc << nl
            << "    loop total      " << wyLoop << nl
            << "    unaccounted     " << (wyLoop - acc) << nl << endl;
    }

    Info<< "End\n" << endl;
    return 0;
}


// ************************************************************************* //
