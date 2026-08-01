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
    interIsoWYFoam

Description
    Two-phase incompressible solver with direction-split Weymouth-Yue
    geometric VOF advection (splitAdvector) and mass-momentum consistent
    per-phase momentum transport inside the sweeps (stage 2 of the
    2026-07 rebuild; Weymouth & Yue JCP 229 2010, Kuhn et al. 2023
    Sec. 3.2 for the consistency coupling).

    Per step: the phase momenta q1 = alpha rho1 u, q2 = (1-alpha)
    rho2 u are split from the one-fluid state, alpha and both momenta
    advance together through the alternating sweeps with the same
    geometric fluxes and the frozen-marker compression terms, and the
    advected velocity is recovered as u* = (q1+q2)/rho(alpha^{n+1}).
    The remaining implicit balance (viscous stresses, gravity,
    pressure via the standard PISO projection with the Sp identity)
    is unchanged from the interIsoConserving family, so all
    differences to that solver family live in the advection block.

    Requirements inherited from splitAdvector: 2D uniform axis-aligned
    static mesh, per-sweep face Courant <= 1/2, cyclic and processor
    boundaries fully supported (a hard requirement of the wave
    cases). nOuterCorrectors 1 and nAlphaSubCycles 1 are enforced.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "geometricAdvector.H"
#include "EulerDdtScheme.H"
#include "localEulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "subCycle.H"
#include "immiscibleIncompressibleTwoPhaseGenericMixture.H"
#include "incompressibleInterPhaseTransportModel.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "CorrectPhi.H"
#include "fvcSmooth.H"
#include "interfaceForces.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Two-phase incompressible solver with Weymouth-Yue split"
        " geometric VOF and consistent per-phase momentum transport"
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "initContinuityErrs.H"
    pimpleControl pimple(mesh);
    #include "createTimeControls.H"
    #include "createFields.H"
    #include "initCorrectPhi.H"
    #include "createUfIfPresent.H"

    if (pimple.nCorrPIMPLE() > 1)
    {
        FatalErrorInFunction
            << "interIsoWYFoam supports nOuterCorrectors 1 only"
            << exit(FatalError);
    }

    Info<< "interIsoWYFoam config:"
        << " reconstruction "
        << mesh.solverDict(alpha1.name())
              .getOrDefault<word>("reconstructionScheme", "plicRDF")
        << ", nCorrectors " << pimple.dict().getOrDefault<label>
           (
               "nCorrectors", 1
           )
        << nl << endl;

    #include "CourantNo.H"
    #include "setInitialDeltaT.H"

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "CourantNo.H"
        #include "alphaCourantNo.H"
        #include "setDeltaT.H"

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            // split the phase momenta from the one-fluid state
            q1 == alpha1*rho1*U;
            q2 == (scalar(1) - alpha1)*rho2*U;

            // consistent geometric advection of alpha, q1, q2
            advector.advect(q1, q2);

            #include "rhofs.H"
            rhoPhi = rho1f*advector.alphaPhi()
                   + rho2f*(phi - advector.alphaPhi());

            alpha2 = 1.0 - alpha1;
            mixture.correct();

            Info<< "Phase-1 volume fraction = "
                << alpha1.weightedAverage(mesh.Vsc()).value()
                << "  Min(" << alpha1.name() << ") = "
                << min(alpha1).value()
                << "  Max(" << alpha1.name() << ") = "
                << max(alpha1).value()
                << endl;

            rho == alpha1*rho1 + alpha2*rho2;

            interfaceForce.correct();

            if (pimple.frozenFlow())
            {
                continue;
            }

            // recovery of the advected one-fluid velocity
            Ustar.primitiveFieldRef() =
                (q1.primitiveField() + q2.primitiveField())
               /rho.primitiveField();
            Ustar.correctBoundaryConditions();

            U = Ustar;
            U.correctBoundaryConditions();

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
