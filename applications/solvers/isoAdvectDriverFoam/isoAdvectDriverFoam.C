/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2026 Asim Onder
     \\/     M anipulation  |
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
    isoAdvectDriverFoam

Description
    Kinematic test driver for the Weymouth-Yue advector (stage 1).
    Advects alpha.water in the prescribed, steady velocity field read
    from 0/U; the flux is computed once by linear interpolation. Each
    step prints the exact mass integral, alpha bounds, and the wisp
    counters, giving the conservation and boundedness record the
    stage-1 gates require. No momentum, no pressure.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "isoAdvection.H"

int main(int argc, char *argv[])
{
    argList::addNote("Weymouth-Yue stage-1 kinematic advection driver");

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"

    Info<< "Reading field alpha.water\n" << endl;
    volScalarField alpha1
    (
        IOobject
        (
            "alpha.water",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    surfaceScalarField phi
    (
        IOobject
        (
            "phi",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvc::flux(U)
    );

    zeroField Su;
    zeroField Sp;
    isoAdvection advector(alpha1, phi, U);

    const scalar mass0 =
        gSum(alpha1.primitiveField()*mesh.V().field());
    Info<< "Initial liquid volume = " << mass0 << " m^3" << nl << endl;

    while (runTime.run())
    {
        ++runTime;

        advector.advect(Sp, Su);

        const scalar mass =
            gSum(alpha1.primitiveField()*mesh.V().field());

        Info<< "t = " << runTime.timeName()
            << "  mass rel. drift = " << (mass - mass0)/mass0
            << "  min/max alpha = "
            << gMin(alpha1.primitiveField()) << " "
            << gMax(alpha1.primitiveField())
            << endl;

        runTime.write();
    }

    Info<< "End\n" << endl;
    return 0;
}

// ************************************************************************* //
