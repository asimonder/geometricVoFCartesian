/*---------------------------------------------------------------------------*\
    v4 consistent-subcycling momentum sub-step.

    Included inside the alpha subCycle loop in alphaEqnSubCycle.H, where
    runTime.deltaT() returns the sub-step dt/m. Advances the working
    phase momentum qWork with THIS sub-step's geometric fluxes
    (advector.alphaPhi() of the sub-step just completed) and refreshes
    the donor velocity uPrev = qWork/rho(alpha) afterwards, so mass and
    momentum stay consistent sweep by sweep and the donor state does not
    freeze over the full step. Motivated by the energy-budget
    measurement of 2026-07-17: the whole numerical sink sits in the
    recovery of the merged momentum, scaling with the mass exchanged per
    sweep times the accumulated donor-velocity mismatch.
\*---------------------------------------------------------------------------*/

{
    const surfaceScalarField& alphaPhiL = advector.alphaPhi();
    const surfaceScalarField alphaPhiG("alphaPhiG", phi - alphaPhiL);

    const scalar limiterTheta =
        pimple.dict().getOrDefault<scalar>("mcvofLimiterTheta", 1.5);
    const bool tracing =
        pimple.dict().getOrDefault<bool>("mcvofTracing", true);
    const scalar halfDt =
        tracing ? 0.5*runTime.deltaTValue() : 0.0;
    const scalar maskTol =
        pimple.dict().getOrDefault<scalar>("mcvofPhaseMaskTol", 1e-2);
    const bool unlimitedSlopes =
        pimple.dict().getOrDefault<bool>("mcvofUnlimitedSlopes", true);

    PtrList<volVectorField> uSlopesL
    (
        mcvofSlopes(mesh, uPrev, limiterTheta,
                    alphaPrev, true, maskTol, unlimitedSlopes)
    );
    PtrList<volVectorField> uSlopesG
    (
        mcvofSlopes(mesh, uPrev, limiterTheta,
                    alphaPrev, false, maskTol, unlimitedSlopes)
    );

    const surfaceVectorField uRecL
    (
        mcvofDonorValue(mesh, alphaPhiL, uPrev, uSlopesL, halfDt)
    );
    const surfaceVectorField uRecG
    (
        mcvofDonorValue(mesh, alphaPhiG, uPrev, uSlopesG, halfDt)
    );

    qWork -= runTime.deltaT()
            *fvc::surfaceIntegrate(rho1*uRecL*alphaPhiL + rho2*uRecG*alphaPhiG);

    // refresh the phase state and working velocity for the next sweep
    alphaPrev = alpha1;
    const volScalarField rhoNow
    (
        rho1*alphaPrev + rho2*(1.0 - alphaPrev)
    );
    uPrev.primitiveFieldRef() =
        qWork.primitiveField()/rhoNow.primitiveField();
    uPrev.correctBoundaryConditions();
}

// ************************************************************************* //
