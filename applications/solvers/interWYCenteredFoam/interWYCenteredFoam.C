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
    interWYCenteredFoam

Description
    Full port of the Basilisk centered + conserving arrangement to the
    cyclic uniform-Cartesian zone class (2D, one empty direction):

      1. volume fraction staggered to half-integer time levels (one
         dt/2 alpha advection at the first step where dt plateaus);
      2. per-phase momenta q1 = alpha rho1 u, q2 = (1-alpha) rho2 u
         advected with alpha through the Weymouth-Yue split sweeps
         (splitAdvector) using the prognostic face flux; velocity
         recovered as u = (q1+q2)/rho(alpha_new);
      3. implicit viscous solve on u + dt g (the previous combined
         pressure-gravity gradient), then g removed again
         (centered.h correction(dt)/correction(-dt) pattern);
      4. prognostic face flux rebuilt from the cell velocity plus the
         face acceleration, ONE exact face projection
         laplacian(dt/rho_f, p_rgh) == div(phi), no Rhie-Chow, no
         ddtCorr, no correctors;
      5. cell velocity corrected by the per-direction two-face AVERAGE
         of the balanced face gradient (approximate projection), not
         by least-squares reconstruction.

    Face density weights use 1/interpolate(rho) (the Basilisk
    alpha_f = 1/rho(face-averaged f)); gravity is the well-balanced
    face form -gh snGrad(rho) sharing its discrete operators with the
    pressure gradient, so the hydrostatic state is machine-exact at
    rest by construction.

    Scope: static 2D uniform axis-aligned mesh, cyclic and wall
    boundaries, laminar. Inherited guards from splitAdvector apply.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "geometricAdvector.H"
#include "splitAdvector.H"
#include "immiscibleIncompressibleTwoPhaseGenericMixture.H"
#include "pimpleControl.H"
#include "emptyFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Basilisk centered-cycle port with Weymouth-Yue consistent"
        " momentum transport on cyclic Cartesian zones"
    );

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "initContinuityErrs.H"

    pimpleControl pimple(mesh);

    #include "createTimeControls.H"

    #include "createFields.H"
    #include "createZonalControls.H"

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

        const dimensionedScalar dtDim("dt", dimTime, runTime.deltaTValue());

        #include "staggerAlpha.H"
        #include "alphaAdvection.H"
        #include "momentumRecovery.H"
        #include "viscousStep.H"
        #include "faceFluxProjection.H"
        #include "cellCorrection.H"

        runTime.write();
        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;
    return 0;
}


// ************************************************************************* //
