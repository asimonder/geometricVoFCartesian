/*---------------------------------------------------------------------------*\
    MCVOF explicit momentum advection.

    Zhang, Popinet, Ling (2020), IJMF 130:103362, Eqs. 12-14 and 18-19,
    with isoAdvector's unsplit geometric liquid volume flux replacing the
    direction-split Weymouth-Yue machinery (unsplit consistency per Kuhn,
    Deskos, Sprague 2023, Sec. 3.2).

    Per step, after the alpha advection:
      liquid volume flux  alphaPhiL = advector.alphaPhi()      [m3/s]
      gas volume flux     alphaPhiG = phi - alphaPhiL          [m3/s]
      F_ml = rho1 * u_donor(rec) * alphaPhiL                  (Eq. 18)
      F_mg = rho2 * u_donor(rec) * alphaPhiG                  (Eq. 19)
      q*   = q^n - dt * div(F_ml + F_mg)                      (Eqs. 12+13)
      u*   = q* / rho^{n+1}                                   (Eq. 14)
    with u_donor reconstructed from U^n by donor-cell values plus minmod
    slopes. A single shared velocity reconstruction serves both phase
    fluxes (one-fluid velocity; the phase split matters only through the
    geometric volume fluxes).

    v1 restrictions enforced here: nAlphaSubCycles == 1.
\*---------------------------------------------------------------------------*/

{
    if (nAlphaSubCycles > 1)
    {
        FatalErrorInFunction
            << "interIsoConservingFoam v1 requires nAlphaSubCycles 1"
            << exit(FatalError);
    }

    if (mesh.changing())
    {
        WarningInFunction
            << "MCVOF momentum advection is untested under mesh change"
            << endl;
    }

    // phase volume fluxes of the just-completed alpha step
    const surfaceScalarField& alphaPhiL = advector.alphaPhi();
    const surfaceScalarField alphaPhiG("alphaPhiG", phi - alphaPhiL);

    // v2 switches (PIMPLE dict)
    const scalar limiterTheta =
        pimple.dict().getOrDefault<scalar>("mcvofLimiterTheta", 1.5);
    const bool tracing =
        pimple.dict().getOrDefault<bool>("mcvofTracing", true);
    const scalar halfDt =
        tracing ? 0.5*runTime.deltaTValue() : 0.0;

    // donor-cell reconstruction of U^n with generalized minmod slopes
    // and BCG-lite half-step tracing
    PtrList<volVectorField> uSlopes
    (
        mcvofSlopes(mesh, U.oldTime(), limiterTheta)
    );

    const surfaceVectorField uRecL
    (
        mcvofDonorValue(mesh, alphaPhiL, U.oldTime(), uSlopes, halfDt)
    );
    const surfaceVectorField uRecG
    (
        mcvofDonorValue(mesh, alphaPhiG, U.oldTime(), uSlopes, halfDt)
    );

    // combined phase-momentum update, Zhang Eqs. 12+13
    volVectorField qStar
    (
        IOobject
        (
            "qStar",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        (rho1*alpha1.oldTime() + rho2*(scalar(1) - alpha1.oldTime()))
       *U.oldTime()
    );

    qStar -= runTime.deltaT()
            *fvc::surfaceIntegrate(rho1*uRecL*alphaPhiL + rho2*uRecG*alphaPhiG);

    // velocity recovery, Zhang Eq. 14; rho already holds the new mixture
    // density (set at the end of alphaEqnSubCycle.H)
    Ustar.primitiveFieldRef() = qStar.primitiveField()/rho.primitiveField();
    Ustar.correctBoundaryConditions();

    U = Ustar;
    U.correctBoundaryConditions();
}

// ************************************************************************* //
