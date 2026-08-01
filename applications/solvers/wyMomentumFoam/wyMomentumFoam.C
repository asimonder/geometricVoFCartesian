/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Copyright (C) 2026 Asim Onder
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM (GPL v3 or later, see
    <http://www.gnu.org/licenses/>).

Application
    wyMomentumFoam

Description
    Stage-2 kinematic gate for the Weymouth-Yue advector: alpha and the
    phase momenta q1 = alpha rho1 u, q2 = (1-alpha) rho2 u advected
    consistently in a prescribed steady velocity field; no pressure, no
    viscosity. For uniform initial velocity the recovered
    u = (q1+q2)/rho must stay uniform to round-off at any density
    ratio (the Galilean consistency identity). Prints per step the
    liquid volume drift, the max relative velocity deviation from the
    initial uniform value, and the momentum sum drift.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "splitAdvector.H"

int main(int argc, char *argv[])
{
    argList::addNote("Weymouth-Yue stage-2 consistency gate driver");

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"

    IOdictionary transportPropertiesEarly
    (
        IOobject
        (
            "transportProperties", runTime.constant(), mesh,
            IOobject::MUST_READ, IOobject::NO_WRITE
        )
    );

    // phases convention (phase subdicts) with flat-entry fallback
    word alphaName("alpha.water");
    if (transportPropertiesEarly.found("phases"))
    {
        const wordList phases
        (
            transportPropertiesEarly.lookup("phases")
        );
        alphaName = "alpha." + phases[0];
    }
    Info<< "Reading field " << alphaName << nl << endl;

    volScalarField alpha1
    (
        IOobject
        (
            alphaName, runTime.timeName(), mesh,
            IOobject::MUST_READ, IOobject::AUTO_WRITE
        ),
        mesh
    );

    volVectorField U
    (
        IOobject
        (
            "U", runTime.timeName(), mesh,
            IOobject::MUST_READ, IOobject::AUTO_WRITE
        ),
        mesh
    );

    surfaceScalarField phi
    (
        IOobject
        (
            "phi", runTime.timeName(), mesh,
            IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE
        ),
        fvc::flux(U)
    );

    const IOdictionary& transportProperties = transportPropertiesEarly;
    dimensionedScalar rho1("rho1", dimDensity, 1000);
    dimensionedScalar rho2("rho2", dimDensity, 1);
    if (transportProperties.found("phases"))
    {
        const wordList phases
        (
            transportProperties.lookup("phases")
        );
        rho1 = dimensionedScalar
        (
            "rho", dimDensity,
            transportProperties.subDict(phases[0])
        );
        rho2 = dimensionedScalar
        (
            "rho", dimDensity,
            transportProperties.subDict(phases[1])
        );
    }
    else
    {
        rho1 = dimensionedScalar
        (
            "rho1", dimDensity, transportProperties
        );
        rho2 = dimensionedScalar
        (
            "rho2", dimDensity, transportProperties
        );
    }

    Info<< "wyMomentumFoam: rho1 " << rho1.value()
        << " rho2 " << rho2.value()
        << " (ratio " << rho1.value()/rho2.value() << ")" << nl << endl;

    volVectorField q1
    (
        IOobject
        (
            "q1", runTime.timeName(), mesh,
            IOobject::NO_READ, IOobject::AUTO_WRITE
        ),
        alpha1*rho1*U
    );
    volVectorField q2
    (
        IOobject
        (
            "q2", runTime.timeName(), mesh,
            IOobject::NO_READ, IOobject::AUTO_WRITE
        ),
        (scalar(1) - alpha1)*rho2*U
    );

    splitAdvector advector(alpha1, phi, U);

    const scalar mass0 = gSum(alpha1.primitiveField()*mesh.V().field());
    const vector mom0 =
        gSum((q1.primitiveField() + q2.primitiveField())*mesh.V().field());
    const vector U0 = U.primitiveField().size() ? U.primitiveField()[0] : vector::zero;
    const scalar magU0 = mag(U0) + VSMALL;

    Info<< "Initial liquid volume " << mass0
        << ", momentum " << mom0 << nl << endl;

    while (runTime.run())
    {
        ++runTime;

        advector.advect(q1, q2);

        // recover the one-fluid velocity
        volScalarField rho(alpha1*rho1 + (scalar(1) - alpha1)*rho2);
        U.primitiveFieldRef() =
            (q1.primitiveField() + q2.primitiveField())
           /rho.primitiveField();
        U.correctBoundaryConditions();

        const scalar mass = gSum(alpha1.primitiveField()*mesh.V().field());
        const vector mom =
            gSum((q1.primitiveField() + q2.primitiveField())*mesh.V().field());

        scalar uDev = 0;
        forAll(U.primitiveField(), celli)
        {
            uDev = max(uDev, mag(U.primitiveField()[celli] - U0));
        }
        reduce(uDev, maxOp<scalar>());

        Info<< "t = " << runTime.timeName()
            << "  mass drift = " << (mass - mass0)/mass0
            << "  max |u-U0|/|U0| = " << uDev/magU0
            << "  mom drift = " << mag(mom - mom0)/(mag(mom0) + VSMALL)
            << "  min/max alpha = " << gMin(alpha1.primitiveField())
            << " " << gMax(alpha1.primitiveField())
            << endl;

        runTime.write();
    }

    Info<< "End\n" << endl;
    return 0;
}

// ************************************************************************* //
