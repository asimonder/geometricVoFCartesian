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
    interIsoConserving4Foam


Description
    Mass-momentum consistent variant (MCVOF), Basilisk-parity pass.
    As interIsoConservingFoam (Zhang, Popinet, Ling 2020, IJMF
    130:103362, Eqs. 12-19) with the v3 per-phase donor reconstruction
    (two slope sets with phase masking, mcvofPhaseMaskTol default 1e-2,
    single-sided at phase edges, shared mixed-cell base value) plus two
    v4 elements, each switchable in the PIMPLE dict:
      mcvofStaggerAlpha (default true): alpha is shifted to the
        half-integer time levels by ONE extra dt/2 advection at a fresh
        start (t = 0 only; restarts inherit the staggering). Every
        subsequent step then advects alpha and momentum over
        [n-1/2, n+1/2] with the projected flux at level n, which sits
        at the midpoint of the advection window. This reproduces the
        time-staggered arrangement of Basilisk (centered.h: tracers lag
        the velocity by dt/2) and the staggered-density formulation of
        Zhang Sec. 2.2.3 WITHOUT any flux extrapolation: the momentum
        continues to use exactly the geometric fluxes isoAdvector
        produces. Diagnosis motivating this: the 2026-07 N512 Lamb
        ladder measured a first-order-in-dt dissipation channel
        r = r0 + b*dt with b = 10.6 1/s^2 attributed to the forward
        (level-n) placement of the advection window.
      mcvofUnlimitedSlopes (default true): donor slopes are plain
        centered differences (the mean of the two one-sided slopes)
        wherever both phase-masked neighbours qualify, instead of
        generalized minmod. This is the reconstruction the Basilisk
        reference actually runs on this benchmark (vof.h
        concentration gradient with u.x.gradient unset = unlimited
        centered in the bulk, one-sided within one cell of the
        interface). Minmod clips to zero slope at every velocity
        extremum (wave crest, orbital maximum), which is measured as
        the first-order spatial dissipation floor (corrected order
        0.975). With mcvofUnlimitedSlopes false the generalized
        minmod path with mcvofLimiterTheta (default 1.5) is retained.
      mcvofTracing (default true): donor values traced back by
        0.5*dt*U_donor (BCG-lite: Taylor term only, no pressure or
        viscous source terms of full BCG - a documented deviation).
    The v2 switch mcvofExtrapolatePhi is REMOVED in v4: the line
    phi = 1.5*phi - 0.5*phi.oldTime() was proven a value-level no-op
    (oldTime() stores the CURRENT field on first access after
    ++runTime; verified in GeometricField.C of v2506), and half-step
    flux extrapolation is superseded by mcvofStaggerAlpha.
    Viscous terms remain implicit Euler and pressure remains at n+1
    (the flat Basilisk reference also runs implicit viscosity).
    With mcvofStaggerAlpha false and mcvofUnlimitedSlopes false, v4
    reproduces v3 exactly (regression control).
    v4: static mesh; nOuterCorrectors 1; nAlphaSubCycles 1.

    Solver derived from interFoam for two incompressible, isothermal immiscible
    fluids using the isoAdvector phase-fraction based interface capturing
    approach, with optional mesh motion and mesh topology changes including
    adaptive re-meshing.

    Reference:
    \verbatim
        Roenby, J., Bredmose, H. and Jasak, H. (2016).
        A computational method for sharp interface advection
        Royal Society Open Science, 3
        doi 10.1098/rsos.160405
    \endverbatim

    isoAdvector code supplied by Johan Roenby, STROMNING (2018)

    Adapted by A. Onder to include new surface tension model (2021-2024)
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
#include "mcvofSlopes.H"
#include "dynamicMultiDimRefineFvMesh.H"
#include "dynamicMultiDimRefineBalancedFvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Solver for two incompressible, isothermal immiscible fluids"
        " using isoAdvector phase-fraction based interface capturing.\n"
        "With optional mesh motion and mesh topology changes including"
        " adaptive re-meshing.\n"
        "The solver is derived from interFoam"
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

    if (pimple.nCorrPIMPLE() > 1)
    {
        FatalErrorInFunction
            << "interIsoConserving4Foam supports nOuterCorrectors 1 only"
            << exit(FatalError);
    }

    if (pimple.dict().found("mcvofExtrapolatePhi"))
    {
        FatalErrorInFunction
            << "mcvofExtrapolatePhi was removed in v4 (proven no-op in "
            << "v1-v3; superseded by mcvofStaggerAlpha)"
            << exit(FatalError);
    }

    Info<< "mcvof config:"
        << " staggerAlpha="
        << pimple.dict().getOrDefault<bool>("mcvofStaggerAlpha", true)
        << " unlimitedSlopes="
        << pimple.dict().getOrDefault<bool>("mcvofUnlimitedSlopes", true)
        << " tracing="
        << pimple.dict().getOrDefault<bool>("mcvofTracing", true)
        << " limiterTheta="
        << pimple.dict().getOrDefault<scalar>("mcvofLimiterTheta", 1.5)
        << " phaseMaskTol="
        << pimple.dict().getOrDefault<scalar>("mcvofPhaseMaskTol", 1e-2)
        << " nAlphaSubCycles="
        << mesh.solverDict(alpha1.name())
              .getOrDefault<label>("nAlphaSubCycles", 1)
        << " (consistent subcycling when > 1)"
        << nl << endl;

    //turbulence->validate();

    #include "CourantNo.H"
    #include "setInitialDeltaT.H"

    // v4: stagger alpha to the half-integer time levels by one extra
    // dt/2 alpha advection, executed ONCE at the first step where the
    // adjustable time step has plateaued (cases start at deltaT 1e-12
    // and ramp up by 1.2x per step; a shift during the ramp is
    // meaningless, which is how the first implementation failed).
    // Advecting with the halved flux over deltaT displaces the
    // interface by deltaT/2; binary halving/doubling restores phi
    // exactly. Fresh starts only; restarts inherit the staggering.
    const bool mcvofStaggerAlpha =
        pimple.dict().getOrDefault<bool>("mcvofStaggerAlpha", true);
    bool mcvofStaggerDone = (runTime.value() > SMALL);
    scalar mcvofDtPrev = -1.0;

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readDyMControls.H"
        #include "CourantNo.H"
        #include "alphaCourantNo.H"
        #include "setDeltaT.H"

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        if (mcvofStaggerAlpha && !mcvofStaggerDone)
        {
            const scalar dtNow = runTime.deltaTValue();
            if (mcvofDtPrev > 0 && mag(dtNow - mcvofDtPrev) < 1e-3*dtNow)
            {
                Info<< "mcvof: dt plateau reached; staggering alpha by"
                    << " dt/2 = " << 0.5*dtNow << " s"
                    << " (one advection with the halved flux)" << nl;

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

                Info<< "mcvof: post-stagger phase-1 volume fraction = "
                    << alpha1.weightedAverage(mesh.Vsc()).value()
                    << "  Min(alpha1) = " << min(alpha1).value()
                    << "  Max(alpha1) = " << max(alpha1).value()
                    << nl << endl;

                mcvofStaggerDone = true;
            }
            mcvofDtPrev = dtNow;
        }

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

            #include "alphaControls.H"

            // v4: no flux manipulation. With mcvofStaggerAlpha the
            // level-n projected flux sits at the midpoint of the
            // [n-1/2, n+1/2] alpha window; the momentum uses exactly
            // the geometric fluxes isoAdvector produces from it.
            // With nAlphaSubCycles > 1 the momentum advances inside
            // the subcycle loop (consistent subcycling) and
            // mcvofAdvection.H below self-skips.
            bool mcvofSubcycled = false;

            #include "alphaEqnSubCycle.H"

            mixture.correct();

	    interfaceForce.correct();

            if (pimple.frozenFlow())
            {
                continue;
            }

            #include "mcvofAdvection.H"

            #include "UEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
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
