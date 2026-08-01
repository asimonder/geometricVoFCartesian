Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh", runTime.timeName(), mesh,
        IOobject::MUST_READ, IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U", runTime.timeName(), mesh,
        IOobject::MUST_READ, IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhi.H"

IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties", runTime.constant(), mesh,
        IOobject::MUST_READ_IF_MODIFIED, IOobject::NO_WRITE
    )
);

Info<< "Reading transportProperties\n" << endl;
immiscibleIncompressibleTwoPhaseGenericMixture mixture(U, phi);

volScalarField& alpha1(mixture.alpha1());
volScalarField& alpha2(mixture.alpha2());

const dimensionedScalar& rho1 = mixture.rho1();
const dimensionedScalar& rho2 = mixture.rho2();

volScalarField rho
(
    IOobject
    (
        "rho", runTime.timeName(), mesh,
        IOobject::READ_IF_PRESENT
    ),
    alpha1*rho1 + alpha2*rho2
);
rho.oldTime();

// air-side driving force (wind forcing) delivered through the
// balanced face acceleration in pEqn; zero vector = inactive
const dimensionedVector airForce
(
    "airForce", dimensionSet(1, -2, -2, 0, 0), vector::zero,
    transportProperties
);
const bool airForcing = mag(airForce.value()) > VSMALL;
if (airForcing)
{
    Info<< "isoCentered: air-side face forcing airForce = "
        << airForce.value() << " [N/m3]" << nl;
}

// turbulence transport: constructed only when turbulenceProperties
// selects a model; laminar cases keep the hand-written viscous path
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
    Info<< "isoCentered: turbulence transport ACTIVE"
        << " (mass-consistent rhoPhi from the geometric alpha flux)"
        << nl;
}

// surface tension: active when sigma > 0; curvature via curvatureDict
const dimensionedScalar sigmaDim
(
    "sigma", dimMass/sqr(dimTime), 0, transportProperties
);
const scalar sigmaVal = sigmaDim.value();
const bool stActive = sigmaVal > VSMALL;
autoPtr<IOdictionary> curvatureDictPtr;
autoPtr<interfaceForces> interfaceForcePtr;
scalar dtCap = GREAT;
if (stActive)
{
    curvatureDictPtr.reset
    (
        new IOdictionary
        (
            IOobject
            (
                "curvatureDict", runTime.system(), mesh,
                IOobject::MUST_READ_IF_MODIFIED, IOobject::NO_WRITE
            )
        )
    );
    interfaceForcePtr.reset
    (
        new interfaceForces(alpha1, U, transportProperties, *curvatureDictPtr)
    );
    // capillary step ceiling: sqrt(rho_mean h^3 / (pi sigma))
    const scalar hMin =
        1.0/gMax(mesh.surfaceInterpolation::deltaCoeffs().primitiveField());
    dtCap = Foam::sqrt
    (
        0.5*(rho1.value() + rho2.value())*Foam::pow3(hMin)
       /(Foam::constant::mathematical::pi*sigmaVal)
    );
    Info<< "interIsoConsistentCenteredFoam: surface tension ACTIVE, sigma = "
        << sigmaVal << ", curvature model "
        << curvatureDictPtr->get<word>("curvatureModel")
        << ", capillary dt cap = " << dtCap << " s" << nl << endl;
}
else
{
    Info<< "interIsoConsistentCenteredFoam: surface tension inactive"
           " (sigma = 0)" << nl << endl;
}

#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"

volScalarField p
(
    IOobject
    (
        "p", runTime.timeName(), mesh,
        IOobject::NO_READ, IOobject::AUTO_WRITE
    ),
    p_rgh + rho*gh
);

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell
(
    p,
    p_rgh,
    pimple.dict(),
    pRefCell,
    pRefValue
);

if (p_rgh.needReference())
{
    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pRefValue - getRefCellValue(p, pRefCell)
    );
    p_rgh = p - rho*gh;
}

mesh.setFluxRequired(p_rgh.name());
mesh.setFluxRequired(alpha1.name());

#include "createFvOptions.H"

isoAdvection advector(alpha1, phi, U);

// MCVOF: advected-state velocity u* (Zhang Eq. 14), source of the
// implicit viscous step in UEqn.H
volVectorField Ustar
(
    IOobject
    (
        "Ustar", runTime.timeName(), mesh,
        IOobject::NO_READ, IOobject::NO_WRITE
    ),
    U
);

// combined pressure-gravity acceleration at cell centres, zero before
// the first projection; persisted for exact restarts
volVectorField gC
(
    IOobject
    (
        "gCentered", runTime.timeName(), mesh,
        IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector(dimAcceleration, Zero)
);

// step-entry velocity for the prognostic-flux proposal (Dev switch)
volVectorField Upre
(
    IOobject
    (
        "Upre", runTime.timeName(), mesh,
        IOobject::NO_READ, IOobject::NO_WRITE, IOobject::NO_REGISTER
    ),
    U
);

// per-direction face-to-cell average (cellCorrection average branch;
// uniform axis-aligned Cartesian meshes only). Mesh accessors are
// fetched per call: topology change reallocates the addressing.
auto faceToCellAvg = [&mesh]
(
    const surfaceScalarField& gf,
    volVectorField& gcell
)
{
    const labelUList& own = mesh.owner();
    const labelUList& nei = mesh.neighbour();
    const surfaceVectorField& Sf = mesh.Sf();
    const surfaceScalarField& magSf = mesh.magSf();

    gcell == dimensionedVector(dimAcceleration, Zero);
    vectorField& gc = gcell.primitiveFieldRef();

    forAll(own, facei)
    {
        const vector nHat = Sf[facei]/magSf[facei];
        for (direction d = 0; d < 3; ++d)
        {
            if (mag(nHat[d]) < 0.99) continue;
            const scalar v = gf[facei]*sign(nHat[d]);
            gc[own[facei]][d] += 0.5*v;
            gc[nei[facei]][d] += 0.5*v;
        }
    }
    forAll(mesh.boundary(), patchi)
    {
        const fvPatch& patch = mesh.boundary()[patchi];
        if (isA<emptyFvPatch>(patch)) continue;
        const labelUList& fc = patch.faceCells();
        const vectorField& pSf = Sf.boundaryField()[patchi];
        const scalarField& pMagSf = magSf.boundaryField()[patchi];
        const scalarField& pGf = gf.boundaryField()[patchi];
        forAll(fc, fi)
        {
            const vector nHat = pSf[fi]/pMagSf[fi];
            for (direction d = 0; d < 3; ++d)
            {
                if (mag(nHat[d]) < 0.99) continue;
                gc[fc[fi]][d] += 0.5*pGf[fi]*sign(nHat[d]);
            }
        }
    }
    gcell.correctBoundaryConditions();
};

const bool staggerAlpha =
    pimple.dict().getOrDefault<bool>("staggerAlpha", true);
bool staggerDone = (runTime.value() > SMALL);
scalar dtPrev = -1.0;

const word cellCorrection
(
    pimple.dict().getOrDefault<word>("cellCorrection", "reconstruct")
);
if (cellCorrection != "reconstruct" && cellCorrection != "average")
{
    FatalErrorInFunction
        << "cellCorrection must be reconstruct or average, got "
        << cellCorrection << exit(FatalError);
}

Info<< "interIsoConsistentCenteredFoam config:"
    << " advector isoAdvector, reconstruction "
    << mesh.solverDict(alpha1.name())
          .getOrDefault<word>("reconstructionScheme", "plicRDF")
    << ", staggerAlpha=" << staggerAlpha
    << ", cellCorrection=" << cellCorrection
    << " (MCVOF consistent momentum, single projection,"
    << " harmonic face-density weights)"
    << nl << endl;
