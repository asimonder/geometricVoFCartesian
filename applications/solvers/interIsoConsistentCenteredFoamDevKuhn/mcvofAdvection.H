/*---------------------------------------------------------------------------*\
    MCVOF explicit momentum advection, DevKuhn variant.

    Base scheme as in interIsoConsistentCenteredFoam (Zhang, Popinet,
    Ling 2020, Eqs. 12-14 and 18-19, with isoAdvector's unsplit
    geometric liquid volume flux). The variant replaces the two
    independent per-phase face velocities by ONE Favre-averaged face
    velocity in the sense of Kuhn, Deskos, Sprague (2023):

      u_f = (rho1 alphaPhiL uRecL + rho2 alphaPhiG uRecG)
           /(rho1 alphaPhiL + rho2 alphaPhiG)

    The denominator is the mass actually crossing the face during the
    step, taken from isoAdvector's own sub-face phase-volume
    integrals, never an independent extrapolation; this removes the
    unbounded-ratio failure of the extrapolated Favre construction
    (forensic record: wetting-front cells at alpha = 0 receiving
    roughly three times the local fluid speed per event). The blend
    weight is clamped to [0,1] and the resulting face velocity is
    bounded componentwise by the two adjacent old-time cell
    velocities; where the transported mass vanishes the donor value
    upwinded by phi is used. Both phase momentum fluxes then share
    this single face velocity:

      F_ml = rho1 uFavre alphaPhiL,  F_mg = rho2 uFavre alphaPhiG.

    Switch: mcvofFavreKuhn (PIMPLE dict, default true). false
    restores the per-phase base scheme for A/B comparison.
\*---------------------------------------------------------------------------*/

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

    // PER-PHASE donor reconstruction: two slope sets with phase
    // masking so the gas flux extrapolates from gas-bearing
    // neighbours and the liquid flux from liquid-bearing ones
    // (single-sided at phase edges). Shared mixed-cell base value.
    const scalar maskTol =
        pimple.dict().getOrDefault<scalar>("mcvofPhaseMaskTol", 1e-2);

    // centered (unlimited) two-sided slopes as in the Basilisk
    // reference; false restores the generalized-minmod path
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

    // ---- Favre face velocity from the flux's own phase partition ----
    const bool favreKuhn =
        pimple.dict().getOrDefault<bool>("mcvofFavreKuhn", true);

    surfaceVectorField uFavre(uRecL);
    if (favreKuhn)
    {
        const labelUList& own = mesh.owner();
        const labelUList& nei = mesh.neighbour();
        const vectorField& Uo = U.oldTime().primitiveField();
        const scalarField& aLf = alphaPhiL.primitiveField();
        const scalarField& aGf = alphaPhiG.primitiveField();
        const scalarField& phif = phi.primitiveField();
        const vectorField& uLf = uRecL.primitiveField();
        const vectorField& uGf = uRecG.primitiveField();
        vectorField& uFf = uFavre.primitiveFieldRef();
        const scalar r1 = rho1.value();
        const scalar r2 = rho2.value();

        forAll(uFf, facei)
        {
            const scalar mL = r1*aLf[facei];
            const scalar mG = r2*aGf[facei];
            const scalar m = mL + mG;
            const vector& ua = Uo[own[facei]];
            const vector& ub = Uo[nei[facei]];
            vector u;
            if (mag(m) > VSMALL)
            {
                scalar w = mL/m;
                w = min(max(w, scalar(0)), scalar(1));
                u = w*uLf[facei] + (scalar(1) - w)*uGf[facei];
            }
            else
            {
                u = (phif[facei] >= 0 ? ua : ub);
            }
            for (direction d = 0; d < 3; ++d)
            {
                const scalar lo = min(ua[d], ub[d]);
                const scalar hi = max(ua[d], ub[d]);
                u[d] = min(max(u[d], lo), hi);
            }
            uFf[facei] = u;
        }

        forAll(mesh.boundary(), patchi)
        {
            const fvPatchVectorField& Uop =
                U.oldTime().boundaryField()[patchi];
            const vectorField Ui(Uop.patchInternalField());
            vectorField Un(Ui);
            if (Uop.coupled())
            {
                Un = refCast<const coupledFvPatchField<vector>>(Uop)
                    .patchNeighbourField();
            }
            const scalarField& paL = alphaPhiL.boundaryField()[patchi];
            const scalarField& paG = alphaPhiG.boundaryField()[patchi];
            const scalarField& pphi = phi.boundaryField()[patchi];
            const vectorField& puL = uRecL.boundaryField()[patchi];
            const vectorField& puG = uRecG.boundaryField()[patchi];
            vectorField& puF = uFavre.boundaryFieldRef()[patchi];

            forAll(puF, fi)
            {
                const scalar mL = r1*paL[fi];
                const scalar mG = r2*paG[fi];
                const scalar m = mL + mG;
                vector u;
                if (mag(m) > VSMALL)
                {
                    scalar w = mL/m;
                    w = min(max(w, scalar(0)), scalar(1));
                    u = w*puL[fi] + (scalar(1) - w)*puG[fi];
                }
                else
                {
                    u = (pphi[fi] >= 0 ? Ui[fi] : Un[fi]);
                }
                for (direction d = 0; d < 3; ++d)
                {
                    const scalar lo = min(Ui[fi][d], Un[fi][d]);
                    const scalar hi = max(Ui[fi][d], Un[fi][d]);
                    u[d] = min(max(u[d], lo), hi);
                }
                puF[fi] = u;
            }
        }
    }

    // per-phase momentum update, carried as two fields so the wisp
    // parts can be zeroed before the density division
    volVectorField q1Star
    (
        IOobject
        (
            "q1Star", runTime.timeName(), mesh,
            IOobject::NO_READ, IOobject::NO_WRITE, IOobject::NO_REGISTER
        ),
        rho1*alpha1.oldTime()*U.oldTime()
    );
    volVectorField q2Star
    (
        IOobject
        (
            "q2Star", runTime.timeName(), mesh,
            IOobject::NO_READ, IOobject::NO_WRITE, IOobject::NO_REGISTER
        ),
        rho2*(scalar(1) - alpha1.oldTime())*U.oldTime()
    );

    if (favreKuhn)
    {
        q1Star -=
            runTime.deltaT()*fvc::surfaceIntegrate(rho1*uFavre*alphaPhiL);
        q2Star -=
            runTime.deltaT()*fvc::surfaceIntegrate(rho2*uFavre*alphaPhiG);
    }
    else
    {
        q1Star -=
            runTime.deltaT()*fvc::surfaceIntegrate(rho1*uRecL*alphaPhiL);
        q2Star -=
            runTime.deltaT()*fvc::surfaceIntegrate(rho2*uRecG*alphaPhiG);
    }

    // wisp cleanup: a phase with content below tolerance in a cell
    // carries no momentum there
    {
        const scalar wispTol =
            pimple.dict().getOrDefault<scalar>("mcvofWispTol", 1e-10);
        const scalarField& aNew = alpha1.primitiveField();
        vectorField& q1f = q1Star.primitiveFieldRef();
        vectorField& q2f = q2Star.primitiveFieldRef();
        forAll(aNew, celli)
        {
            if (aNew[celli] < wispTol)
            {
                q1f[celli] = vector::zero;
            }
            if (aNew[celli] > scalar(1) - wispTol)
            {
                q2f[celli] = vector::zero;
            }
        }
    }

    // velocity recovery; rho already holds the new mixture density
    Ustar.primitiveFieldRef() =
        (q1Star.primitiveField() + q2Star.primitiveField())
       /rho.primitiveField();
    Ustar.correctBoundaryConditions();

    U = Ustar;
    U.correctBoundaryConditions();
}

// ************************************************************************* //
