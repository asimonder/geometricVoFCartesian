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
    with u_donor reconstructed from U^n by donor-cell values plus
    per-phase directional slopes (centered when mcvofUnlimitedSlopes,
    else generalized minmod). With mcvofStaggerAlpha the alpha levels
    are the half-integer times, so alpha1.oldTime() below is
    alpha^{n-1/2} and the update reads
      q*  = rho(alpha^{n-1/2}) u^n - dt div(F)     (Basilisk pairing)
      u*  = q* / rho(alpha^{n+1/2})
    which matches conserving.h (q1 = f rho1 u with f lagging u by
    dt/2 and recovery through rho(f^{new})).

    v1 restrictions enforced here: nAlphaSubCycles == 1.
\*---------------------------------------------------------------------------*/

if (!mcvofSubcycled)
{
    if (mesh.changing())
    {
        WarningInFunction
            << "MCVOF momentum advection is untested under mesh change"
            << endl;
    }

    // phase volume fluxes of the just-completed alpha step
    const surfaceScalarField& alphaPhiL = advector.alphaPhi();
    const surfaceScalarField alphaPhiG("alphaPhiG", phi - alphaPhiL);

    // switches (PIMPLE dict)
    const scalar limiterTheta =
        pimple.dict().getOrDefault<scalar>("mcvofLimiterTheta", 1.5);
    const bool tracing =
        pimple.dict().getOrDefault<bool>("mcvofTracing", true);
    const scalar halfDt =
        tracing ? 0.5*runTime.deltaTValue() : 0.0;

    // v3: PER-PHASE donor reconstruction. Two slope sets with phase
    // masking so the gas flux extrapolates from gas-bearing neighbours
    // and the liquid flux from liquid-bearing ones (single-sided at
    // phase edges). Shared mixed-cell base value as in the reference.
    const scalar maskTol =
        pimple.dict().getOrDefault<scalar>("mcvofPhaseMaskTol", 1e-2);

    // v4: centered (unlimited) two-sided slopes as in the Basilisk
    // reference; false restores the v3 generalized-minmod path.
    const bool unlimitedSlopes =
        pimple.dict().getOrDefault<bool>("mcvofUnlimitedSlopes", true);

    PtrList<volVectorField> uSlopesL
    (
        mcvofSlopes(mesh, U.oldTime(), limiterTheta,
                    alpha1.oldTime(), true, maskTol, unlimitedSlopes)
    );
    PtrList<volVectorField> uSlopesG
    (
        mcvofSlopes(mesh, U.oldTime(), limiterTheta,
                    alpha1.oldTime(), false, maskTol, unlimitedSlopes)
    );

    const surfaceVectorField uRecL
    (
        mcvofDonorValue(mesh, alphaPhiL, U.oldTime(), uSlopesL, halfDt)
    );
    const surfaceVectorField uRecG
    (
        mcvofDonorValue(mesh, alphaPhiG, U.oldTime(), uSlopesG, halfDt)
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
