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
    (scalar(1) - alpha1)*rho2*U
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

// One-sided cell correction beside the interface.  Per direction, a face
// touching a mixed cell is dropped from the cell's stencil and the
// remaining face carries the whole weight.  Identical to the parent
// wherever a direction's two faces are both clean or both touch a mixed
// cell.  Default false = parent path.
const bool cellCorrectionOneSided
(
    pimple.dict().getOrDefault<Switch>("cellCorrectionOneSided", false)
);
const scalar cellCorrectionAlphaTol
(
    pimple.dict().getOrDefault<scalar>("cellCorrectionAlphaTol", 1e-6)
);
// 1 = adjacent clean face at full weight (first order)
// 2 = one-sided second-order difference from the two clean faces
const label cellCorrectionOneSidedOrder
(
    pimple.dict().getOrDefault<label>("cellCorrectionOneSidedOrder", 1)
);
if (cellCorrectionOneSided)
{
    Info<< "jump: one-sided cell correction beside the interface"
        << " (cellCorrectionOneSided true, alphaTol "
        << cellCorrectionAlphaTol << ")" << nl;
}

// per-direction face-to-cell average that drops the faces touching a
// mixed cell.  Same 0.5 weights as the parent, rescaled by the number of
// faces actually used, so a direction with two clean faces is identical.
auto faceToCellAvgOneSided =
[&mesh, &alpha1, &cellCorrectionAlphaTol, &cellCorrectionOneSidedOrder]
(
    const surfaceScalarField& gfld,
    volVectorField& gcell
)
{
    const labelUList& own = mesh.owner();
    const labelUList& nei = mesh.neighbour();
    const surfaceVectorField& Sf = mesh.Sf();
    const surfaceScalarField& magSf = mesh.magSf();
    const scalarField& a = alpha1.primitiveField();
    const scalar aTol = cellCorrectionAlphaTol;

    gcell == dimensionedVector(dimAcceleration, Zero);
    vectorField& gc = gcell.primitiveFieldRef();
    vectorField cleanSum(gc.size(), Zero);
    vectorField nAll(gc.size(), Zero);
    vectorField nClean(gc.size(), Zero);
    // cell across the clean face of each direction, for the second-order
    // one-sided form; only read where a direction has exactly one clean
    // face, so a cell with two clean faces overwriting it is harmless
    labelList cleanNbr(3*gc.size(), -1);

    auto isMixed = [&a, aTol](const label c)
    {
        return a[c] > aTol && a[c] < 1 - aTol;
    };
    auto add =
    [&gc, &cleanSum, &nAll, &nClean]
    (const label c, const direction d, const scalar v, const bool clean)
    {
        gc[c][d] += 0.5*v;
        nAll[c][d] += 1;
        if (clean)
        {
            cleanSum[c][d] += 0.5*v;
            nClean[c][d] += 1;
        }
    };

    forAll(own, facei)
    {
        const vector nHat = Sf[facei]/magSf[facei];
        const bool clean = !isMixed(own[facei]) && !isMixed(nei[facei]);
        for (direction d = 0; d < 3; ++d)
        {
            if (mag(nHat[d]) < 0.99) continue;
            const scalar v = gfld[facei]*sign(nHat[d]);
            add(own[facei], d, v, clean);
            add(nei[facei], d, v, clean);
            if (clean)
            {
                cleanNbr[3*own[facei] + d] = nei[facei];
                cleanNbr[3*nei[facei] + d] = own[facei];
            }
        }
    }

    forAll(mesh.boundary(), patchi)
    {
        const fvPatch& patch = mesh.boundary()[patchi];
        if (isA<emptyFvPatch>(patch)) continue;

        const labelUList& fc = patch.faceCells();
        const vectorField& pSf = Sf.boundaryField()[patchi];
        const scalarField& pMagSf = magSf.boundaryField()[patchi];
        const scalarField& pGf = gfld.boundaryField()[patchi];

        tmp<scalarField> taNei;
        if (patch.coupled())
        {
            taNei = alpha1.boundaryField()[patchi].patchNeighbourField();
        }

        forAll(fc, fi)
        {
            bool clean = !isMixed(fc[fi]);
            if (clean && patch.coupled())
            {
                const scalar an = taNei()[fi];
                clean = !(an > aTol && an < 1 - aTol);
            }
            const vector nHat = pSf[fi]/pMagSf[fi];
            for (direction d = 0; d < 3; ++d)
            {
                if (mag(nHat[d]) < 0.99) continue;
                add(fc[fi], d, pGf[fi]*sign(nHat[d]), clean);
            }
        }
    }

    forAll(gc, celli)
    {
        for (direction d = 0; d < 3; ++d)
        {
            const scalar na = nAll[celli][d];
            const scalar nc = nClean[celli][d];
            if (!(na > 0 && nc > 0 && nc < na)) continue;

            // v1, the adjacent clean face value at full weight
            const scalar v1 = cleanSum[celli][d]*(na/nc);
            gc[celli][d] = v1;

            if (cellCorrectionOneSidedOrder < 2) continue;

            const label nb = cleanNbr[3*celli + d];
            if (nb < 0 || nClean[nb][d] < 1.5) continue;

            // the neighbour holds S = (v1 + v2)/2 with both faces clean
            gc[celli][d] = 2.0*v1 - cleanSum[nb][d];
        }
    }
    gcell.correctBoundaryConditions();
};

const bool staggerAlpha =
    pimple.dict().getOrDefault<bool>("wyStaggerAlpha", true);

// Interfacial jump treatment of the face density in the pressure equation.
// At a face the interface crosses, the 1D-exact Poisson coefficient is the
// series resistance of the two phase segments,
//
//     1/rho_f = lambda/rho_s(P) + (1 - lambda)/rho_s(N),
//
// with lambda the interface position along the owner-neighbour segment and
// s(.) the phase of a cell.  The harmonic mixture average that the parent
// solver uses is this expression with lambda pinned to 1/2 and the mixture
// cell densities substituted for the phase densities
// (Vukcevic et al. 2017, Eqs. 39-49).  Default false = parent path.
const bool jumpFaceDensity
(
    pimple.dict().getOrDefault<Switch>("jumpFaceDensity", false)
);

// re-reconstruct the interface after the advection, so the distance
// function and the density field describe the same interface
const bool jumpReconstruct
(
    pimple.dict().getOrDefault<Switch>("jumpReconstruct", true)
);

const bool jumpReport
(
    pimple.dict().getOrDefault<Switch>("jumpReport", true)
);
const label jumpReportInterval
(
    pimple.dict().getOrDefault<label>("jumpReportInterval", 100)
);

// Face density from the density the sweeps advect.  The geometric alpha
// flux defines a face volume fraction alpha_f = phi_alpha/phi, and the
// advected face density is rho1 alpha_f + rho2 (1 - alpha_f).  This is the
// quantity the consistency argument puts in the Poisson coefficient
// (Zuzio et al. 2020, Eqs. 24-26; Kuhn et al. 2023, Eq. 30), in place of
// the parent's harmonic average of the mixture cell densities.
// Default false = parent path.
const bool advectedFaceDensity
(
    pimple.dict().getOrDefault<Switch>("advectedFaceDensity", false)
);

// blend width where phi passes through zero, as a fraction of max|phi|
const scalar advectedPhiTol
(
    pimple.dict().getOrDefault<scalar>("advectedPhiTol", 1e-3)
);

if (advectedFaceDensity && jumpFaceDensity)
{
    FatalErrorInFunction
        << "advectedFaceDensity and jumpFaceDensity are alternative face"
        << " densities; enable at most one" << exit(FatalError);
}
if (advectedFaceDensity)
{
    Info<< "jump: face density from the advected alpha flux"
        << " (advectedFaceDensity true, advectedPhiTol "
        << advectedPhiTol << ")" << nl;
}

// Cell correction from the face acceleration rather than the face force.
// The parent averages hf = rho_f gf and divides by the mixture cell
// density.  In a mixed cell that mixes two face forces separated by the
// density ratio and divides by a density belonging to neither phase.
// The face acceleration gf is the material acceleration and is continuous
// across the interface.  Identical to the parent wherever the density is
// constant.  Default false = parent path.
const bool cellCorrectionAcceleration
(
    pimple.dict().getOrDefault<Switch>("cellCorrectionAcceleration", false)
);
if (cellCorrectionAcceleration)
{
    Info<< "jump: cell correction averages the face acceleration"
        << " (cellCorrectionAcceleration true)" << nl;
}

// Buoyancy distributed by a regularised Heaviside of the distance
// function rather than by the volume-fraction differences.  The total and
// the first moment per column are unchanged; what changes is that the
// distribution now varies continuously as the interface crosses a cell
// instead of jumping.  Sussman, Smereka and Osher (1994) smeared
// Heaviside.  Default false = parent path.
const bool buoyancyKernel
(
    pimple.dict().getOrDefault<Switch>("buoyancyKernel", false)
);
const scalar buoyancyKernelWidth
(
    pimple.dict().getOrDefault<scalar>("buoyancyKernelWidth", 1.0)
);
if (buoyancyKernel)
{
    Info<< "jump: buoyancy distributed by a regularised Heaviside, width "
        << buoyancyKernelWidth << " cells" << nl;
}

// which face density is used at the faces the interface crosses.
//   oneSided : lambda/rho_s(P) + (1-lambda)/rho_s(N), the 1D-exact series
//              value (Vukcevic et al. 2017, Eqs. 39-49)
//   queutey  : the arithmetic mean of the two side densities, rho_hat of
//              Queutey and Visonneau (2007), Eqs. (36b), (37), (38)
const word jumpFaceDensityModel
(
    pimple.dict().getOrDefault<word>("jumpFaceDensityModel", "oneSided")
);
if (jumpFaceDensity)
{
    if (jumpFaceDensityModel != "oneSided"
     && jumpFaceDensityModel != "queutey")
    {
        FatalErrorInFunction
            << "jumpFaceDensityModel must be oneSided or queutey, got "
            << jumpFaceDensityModel << exit(FatalError);
    }
    Info<< "jump: interface face density model " << jumpFaceDensityModel
        << nl;
}

autoPtr<reconstructedDistanceFunction> rdfPtr;
if (jumpFaceDensity || buoyancyKernel)
{
    rdfPtr.reset(new reconstructedDistanceFunction(mesh));
    Info<< "jump: one-sided face density at interface faces"
        << " (jumpFaceDensity true, jumpReconstruct "
        << Switch(jumpReconstruct) << ")" << nl;
}
else
{
    Info<< "jump: face density is the harmonic mixture average"
        << " (jumpFaceDensity false, parent path)" << nl;
}
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
