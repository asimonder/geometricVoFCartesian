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
    interWYConsistentCenteredFoam

Description
    Mass-momentum consistent two-phase solver on collocated grids.
    Structure follows interIsoCartFoam. The pressure-velocity cycle is a
    single exact face projection with an averaged cell correction:

      alphaEqn:  per-phase momenta q1 = alpha rho1 u, q2 = (1-alpha) rho2 u
                 advected with alpha by the Weymouth-Yue split sweeps
                 (splitAdvector); volume fraction staggered dt/2 once at
                 the dt plateau.
      UEqn:      velocity recovered as u = (q1+q2)/rho (consistent
                 transport); implicit viscous step.
      pEqn:      face flux from cell velocity plus the balanced face
                 forces (surface tension, gravity); one exact projection
                 with HARMONIC face density rAUf = dt/rho_f,
                 rho_f = 1/interp(1/rho); cell velocities corrected by
                 the per-direction two-face average of the face force
                 density divided by the CELL density.

    Surface tension uses interfaceForces (curvature via curvatureDict),
    active when sigma > 0; sigma = 0 skips construction entirely.

    Dynamic mesh: mesh.update() hooks as in interIsoCartFoam; the
    Weymouth-Yue zone machinery is not yet validated under topology
    change and warns if the mesh changes.

    nOuterCorrectors and nAlphaSubCycles other than 1 are rejected: the
    consistency of the cycle is structural and admits no correctors.
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "geometricAdvector.H"
#include "splitAdvector.H"
#include "immiscibleIncompressibleTwoPhaseGenericMixture.H"
#include "pimpleControl.H"
#include "interfaceForces.H"
#include "reconstructionSchemes.H"
#include "cutFacePLIC.H"
#include "CorrectPhi.H"
#include "adjustPhi.H"
#include "fixedValueFvPatchFields.H"
#include "mathematicalConstants.H"
#include "emptyFvPatch.H"
#include "fvOptions.H"
#include "incompressibleInterPhaseTransportModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Mass-momentum consistent two-phase solver on collocated grids"
        " (Weymouth-Yue consistent transport, single-projection cycle,"
        " harmonic face-density weights)."
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
