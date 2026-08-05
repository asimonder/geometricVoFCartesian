/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2016 DHI
    Copyright (C) 2017 OpenCFD Ltd.
    Copyright (C) 2018 Johan Roenby
    Copyright (C) 2019-2020 DLR
    Copyright (C) 2021 Asim Onder
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
    interIsoDynCartFoamTimed

Description
    Timed clone of interIsoDynCartFoam. The algorithm is unchanged. The
    only additions are wall-clock timers around the stage includes of the
    time loop, and the report they print.

    The clock is a single clockTime instance. Consecutive calls to
    timeIncrement() partition the step exactly, so the stage totals sum
    to the loop total up to the cost of the calls themselves. The stage
    names match the interWYConsistentCenteredFoamTimed build, so the two
    reports compare line for line.

      ctl    controls, Courant numbers, setDeltaT
      mesh   mesh.update, alpha remap, MRF, correctPhi
      adv    alphaEqnSubCycle.H          (isoAdvector)
      surf   mixture.correct, curvature and surface tension
      visc   UEqn.H
      pres   the whole pimple.correct() loop over pEqn.H
      turb   turbulence->correct
      write  runTime.write

    pres covers every pressure corrector of the step, and npc reports how
    many ran. Divide pres by npc for the cost of one corrector.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "isoAdvection.H"
#include "EulerDdtScheme.H"
#include "localEulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "subCycle.H"
#include "immiscibleIncompressibleTwoPhaseGenericMixture.H"
#include "incompressibleInterPhaseTransportModel.H"
//#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "CorrectPhi.H"
#include "fvcSmooth.H"
#include "dynamicRefineFvMesh.H"
#include "interfaceForces.H"
#include "dynamicMultiDimRefineFvMesh.H"
#include "dynamicMultiDimRefineBalancedFvMesh.H"
#include "clockTime.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Solver for two incompressible, isothermal immiscible fluids"
        " using isoAdvector phase-fraction based interface capturing.\n"
        "With optional mesh motion and mesh topology changes including"
        " adaptive re-meshing.\n"
        "The solver is derived from interFoam. Timed build."
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "initContinuityErrs.H"
    #include "createDyMControls.H"
    #include "createFields.H"
    #include "initCorrectPhi.H"
    #include "createUfIfPresent.H"

    //turbulence->validate();

    #include "CourantNo.H"
    #include "setInitialDeltaT.H"

    // --- stage timers (timed build only) -----------------------------------
    clockTime wyClk;
    clockTime wyLoopClk;
    scalar wyCtl = 0, wyMesh = 0, wyAdv = 0, wySurf = 0,
           wyVisc = 0, wyPres = 0, wyTurb = 0, wyWrite = 0;
    label wySteps = 0, wyPresCorr = 0;
    const label wyReport =
        runTime.controlDict().getOrDefault<label>("wyTimerInterval", 1);
    Info<< "wyTimers: stage timing active, report interval "
        << wyReport << " step(s)" << endl;
    // -----------------------------------------------------------------------

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    Info<< "\nStarting time loop\n" << endl;

    wyLoopClk.timeIncrement();

    while (runTime.run())
    {
        wyClk.timeIncrement();

        #include "readDyMControls.H"
        #include "CourantNo.H"
        #include "alphaCourantNo.H"
        #include "setDeltaT.H"

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        const scalar sCtl = wyClk.timeIncrement();
        scalar sMesh = 0, sAdv = 0, sSurf = 0, sVisc = 0, sPres = 0, sTurb = 0;
        label nPresCorr = 0;

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
	  if (pimple.firstIter() || moveMeshOuterCorrectors)
            {
	      if (isA<dynamicRefineFvMesh>(mesh) || isA<dynamicMultiDimRefineFvMesh>(mesh) || isA<dynamicMultiDimRefineBalancedFvMesh>(mesh))
		{
		  advector.surf().reconstruct();
		}

                mesh.update();

                if (mesh.changing())
                {
		  gh = (g & mesh.C()) - ghRef;
		  ghf = (g & mesh.Cf()) - ghRef;

		  //  if (isA<dynamicRefineFvMesh>(mesh))
		  if (isA<dynamicRefineFvMesh>(mesh) || isA<dynamicMultiDimRefineFvMesh>(mesh) || isA<dynamicMultiDimRefineBalancedFvMesh>(mesh))
                    {
		      advector.surf().mapAlphaField();
		      alpha2 = 1.0 - alpha1;
		      alpha2.correctBoundaryConditions();
		      rho == alpha1*rho1 + alpha2*rho2;
		      rho.correctBoundaryConditions();
		      rho.oldTime() = rho;
		      alpha2.oldTime() = alpha2;
                    }

		  MRF.update();

		  if (correctPhi)
                    {
                        // Calculate absolute flux
                        // from the mapped surface velocity
		      phi = mesh.Sf() & Uf();

                      #include "correctPhi.H"

                        // Make the flux relative to the mesh motion
		      fvc::makeRelative(phi, U);

		      mixture.correct();
                    }

                    if (checkMeshCourantNo)
                    {
                        #include "meshCourantNo.H"
		    }

                }
	}

            sMesh += wyClk.timeIncrement();

            #include "alphaControls.H"
            #include "alphaEqnSubCycle.H"

            sAdv += wyClk.timeIncrement();

            mixture.correct();

	    interfaceForce.correct();

            sSurf += wyClk.timeIncrement();

            if (pimple.frozenFlow())
            {
                continue;
            }

            #include "UEqn.H"

            sVisc += wyClk.timeIncrement();

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
                ++nPresCorr;
            }

            sPres += wyClk.timeIncrement();

            if (pimple.turbCorr())
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
        wyPresCorr += nPresCorr;
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
                << "  npc " << nPresCorr
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
            << "    pressure correctors " << wyPresCorr << nl
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
