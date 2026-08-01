/*---------------------------------------------------------------------------*\
    Part of interWYCenteredFoam (see the main file for the full
    description and license). Extracted for the conventional
    fragmented solver layout.
\*---------------------------------------------------------------------------*/

    // --- fields -------------------------------------------------------

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

    #include "readGravitationalAcceleration.H"
    #include "readhRef.H"
    #include "gh.H"

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
        Info<< "interWYCenteredFoam: zonal mode, "
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

    // combined pressure-gravity acceleration at cell centres
    // (centered.h field g), zero before the first projection
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

