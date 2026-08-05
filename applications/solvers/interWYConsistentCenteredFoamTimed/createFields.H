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
    Info<< "wy: air-side face forcing airForce = " << airForce.value()
        << " [N/m3]" << nl;
}

// turbulence transport (interIsoCartFoam pattern): constructed only
// when turbulenceProperties selects a model; laminar cases keep the
// original viscous path bit for bit (laminar bypass)
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
    Info<< "wy: turbulence transport ACTIVE"
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
    Info<< "interWYConsistentCenteredFoam: surface tension ACTIVE, sigma = "
        << sigmaVal << ", curvature model "
        << curvatureDictPtr->get<word>("curvatureModel")
        << ", capillary dt cap = " << dtCap << " s" << nl << endl;
}
else
{
    Info<< "interWYConsistentCenteredFoam: surface tension inactive"
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
{
    // accept the solver family's pRefPoint/pRefValue convention
    dictionary refDict(pimple.dict());
    if (pimple.dict().found("pRefPoint"))
    {
        refDict.add
        (
            "p_rghRefPoint",
            pimple.dict().get<point>("pRefPoint"),
            true
        );
    }
    if (pimple.dict().found("pRefCell"))
    {
        refDict.add
        (
            "p_rghRefCell",
            pimple.dict().get<label>("pRefCell"),
            true
        );
    }
    refDict.add
    (
        "p_rghRefValue",
        pimple.dict().getOrDefault<scalar>("pRefValue", 0),
        true
    );
    setRefCell(p_rgh, refDict, pRefCell, pRefValue);
}
mesh.setFluxRequired(p_rgh.name());
mesh.setFluxRequired(alpha1.name());

autoPtr<geometricAdvector> advectorPtr
(
    geometricAdvector::New(alpha1, phi, U)
);
geometricAdvector& advector = *advectorPtr;

// zonal operation: outside-zone momentum runs through a standard
// flux-form convection operator stitched at the zone faces
const splitAdvector* splitPtr = isA<splitAdvector>(advector)
    ? &refCast<const splitAdvector>(advector) : nullptr;
const bool zonal = (splitPtr && !splitPtr->wholeDomain());
boolList zoneFace;
if (zonal)
{
    const boolList& zc = splitPtr->inZone();
    zoneFace.setSize(mesh.nInternalFaces(), false);
    forAll(zoneFace, facei)
    {
        zoneFace[facei] =
            zc[mesh.owner()[facei]] || zc[mesh.neighbour()[facei]];
    }
    label nOut = 0;
    forAll(zc, c) if (!zc[c]) ++nOut;
    Info<< "interWYConsistentCenteredFoam: zonal mode, "
        << returnReduce(nOut, sumOp<label>())
        << " cells outside the lattice; momentum there is advected "
        << "by the per-phase Lax-Wendroff sweep fluxes"
        << nl << endl;
}

volVectorField q1
(
    IOobject
    (
        "q1", runTime.timeName(), mesh,
        IOobject::NO_READ, IOobject::NO_WRITE
    ),
    alpha1*rho1*U
);
volVectorField q2
(
    IOobject
    (
        "q2", runTime.timeName(), mesh,
        IOobject::NO_READ, IOobject::NO_WRITE
    ),
    (1.0 - alpha1)*rho2*U
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

// per-direction face-to-cell average (the approximate projection's
// cell correction; uniform Cartesian). Mesh accessors are fetched per
// call: topology change reallocates owner/Sf and captured references
// would dangle under AMR.
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
    pimple.dict().getOrDefault<bool>("wyStaggerAlpha", true);
bool staggerDone = (runTime.value() > SMALL);
scalar dtPrev = -1.0;

Info<< "interWYConsistentCenteredFoam config:"
    << " reconstruction "
    << mesh.solverDict(alpha1.name())
          .getOrDefault<word>("reconstructionScheme", "plicRDF")
    << " staggerAlpha=" << staggerAlpha
    << " (single projection, averaged cell correction,"
    << " harmonic face-density weights)"
    << nl << endl;

#include "createFvOptions.H"
