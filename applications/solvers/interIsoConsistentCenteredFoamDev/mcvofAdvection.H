/*---------------------------------------------------------------------------*\
    MCVOF explicit momentum advection, DevCentroid variant.

    DevKuhn base (bounded Favre face velocity on isoAdvector's own
    flux partition) plus the flux time-centroid pairing: the donor
    velocity of each phase is evaluated at the time centroid of that
    phase's face flux instead of the uniform step midpoint. On a face
    the interface crosses during the step, the phase flux is
    concentrated in a sub-interval, and sampling the velocity at dt/2
    pairs mass with motion at the wrong time; this mismatch is the
    injection channel of the centered arrangement at interface faces.

    The centroid is modeled from a linear wetting history of the
    donor cell. With a0 and a1 the donor volume fractions before and
    after the advection, the centroid fraction of the liquid flux is

        s_L = (a0 + 2 a1) / (3 (a0 + a1)),   clamped to [1/3, 2/3],

    which is exactly 1/2 on steady and pure-phase faces (a0 = a1),
    2/3 on wetting faces, 1/3 on drying faces; the gas value uses the
    complements. The evaluation time handed to the traceback is then
    s*dt per face and phase; pure-phase faces reproduce the DevKuhn
    scheme bit for bit.

    Switches (PIMPLE dict):
      mcvofFluxCentroid  (default true)  per-face centroid times
      mcvofFavreKuhn     (default true)  bounded Favre face velocity
      mcvofTracing       (default true)  traceback active at all

    The traceback is always centered at the step midpoint. An earlier
    mcvofTimeTheta switch scaled the traceback time globally; it was
    removed because it is a spatial displacement of the sampling point
    along the local velocity, not a temporal off-centering, and above
    theta*Co = 1/2 it drives the reconstruction coefficient negative
    and unbounded in the light phase.
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
    const bool fluxCentroid =
        tracing
     && pimple.dict().getOrDefault<bool>("mcvofFluxCentroid", true);
    const scalar halfDt =
        tracing ? 0.5*runTime.deltaTValue() : 0.0;

    // ---- per-face, per-phase flux time-centroid evaluation times ----
    autoPtr<surfaceScalarField> tEvalLPtr;
    autoPtr<surfaceScalarField> tEvalGPtr;
    if (fluxCentroid)
    {
        const scalar dtVal = runTime.deltaTValue();

        tEvalLPtr.reset
        (
            new surfaceScalarField
            (
                IOobject
                (
                    "tEvalL", runTime.timeName(), mesh,
                    IOobject::NO_READ, IOobject::NO_WRITE,
                    IOobject::NO_REGISTER
                ),
                mesh,
                dimensionedScalar(dimTime, 0.5*dtVal)
            )
        );
        tEvalGPtr.reset
        (
            new surfaceScalarField
            (
                IOobject
                (
                    "tEvalG", runTime.timeName(), mesh,
                    IOobject::NO_READ, IOobject::NO_WRITE,
                    IOobject::NO_REGISTER
                ),
                mesh,
                dimensionedScalar(dimTime, 0.5*dtVal)
            )
        );

        const scalar sLo = 1.0/3.0;
        const scalar sHi = 2.0/3.0;

        auto centroidFrac = [sLo,sHi](const scalar b0, const scalar b1)
        {
            const scalar den = b0 + b1;
            return
                (den < SMALL)
              ? scalar(0.5)
              : min(max((b0 + 2.0*b1)/(3.0*den), sLo), sHi);
        };

        const labelUList& own = mesh.owner();
        const labelUList& nei = mesh.neighbour();
        const scalarField& a0 = alpha1.oldTime().primitiveField();
        const scalarField& a1f = alpha1.primitiveField();
        const scalarField& fL = alphaPhiL.primitiveField();
        const scalarField& fG = alphaPhiG.primitiveField();
        scalarField& tL = tEvalLPtr().primitiveFieldRef();
        scalarField& tG = tEvalGPtr().primitiveFieldRef();

        forAll(tL, facei)
        {
            const label dcL = (fL[facei] >= 0) ? own[facei] : nei[facei];
            tL[facei] = dtVal*centroidFrac(a0[dcL], a1f[dcL]);

            const label dcG = (fG[facei] >= 0) ? own[facei] : nei[facei];
            tG[facei] =
                dtVal*centroidFrac
                (
                    scalar(1) - a0[dcG], scalar(1) - a1f[dcG]
                );
        }

        // coupled patches (cyclic and processor seams): remote donor
        // fractions from the neighbour fields; physical patches keep
        // the midpoint value
        forAll(mesh.boundary(), patchi)
        {
            if (!mesh.boundary()[patchi].coupled()) continue;

            const scalarField pA0In
            (
                alpha1.oldTime().boundaryField()[patchi]
                    .patchInternalField()
            );
            const scalarField pA0Nei
            (
                refCast<const coupledFvPatchField<scalar>>
                (
                    alpha1.oldTime().boundaryField()[patchi]
                ).patchNeighbourField()
            );
            const scalarField pA1In
            (
                alpha1.boundaryField()[patchi].patchInternalField()
            );
            const scalarField pA1Nei
            (
                refCast<const coupledFvPatchField<scalar>>
                (
                    alpha1.boundaryField()[patchi]
                ).patchNeighbourField()
            );
            const scalarField& pfL = alphaPhiL.boundaryField()[patchi];
            const scalarField& pfG = alphaPhiG.boundaryField()[patchi];
            scalarField& ptL = tEvalLPtr().boundaryFieldRef()[patchi];
            scalarField& ptG = tEvalGPtr().boundaryFieldRef()[patchi];

            forAll(ptL, fi)
            {
                const bool ownL = (pfL[fi] >= 0);
                ptL[fi] = dtVal*centroidFrac
                (
                    ownL ? pA0In[fi] : pA0Nei[fi],
                    ownL ? pA1In[fi] : pA1Nei[fi]
                );

                const bool ownG = (pfG[fi] >= 0);
                ptG[fi] = dtVal*centroidFrac
                (
                    scalar(1) - (ownG ? pA0In[fi] : pA0Nei[fi]),
                    scalar(1) - (ownG ? pA1In[fi] : pA1Nei[fi])
                );
            }
        }
    }

    // PER-PHASE donor reconstruction: two slope sets with phase
    // masking so the gas flux extrapolates from gas-bearing
    // neighbours and the liquid flux from liquid-bearing ones
    // (single-sided at phase edges). Shared mixed-cell base value.
    const scalar maskTol =
        pimple.dict().getOrDefault<scalar>("mcvofPhaseMaskTol", 0.5);

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

    // Basilisk donor-majority rule (vof_concentration_gradient):
    // a cell whose own content of a phase is minority carries zero
    // slope for that phase, exactly as in the WY port's concGradient.
    // Combine with mcvofPhaseMaskTol 0.5 for the faithful Basilisk
    // reconstruction. Default off preserves current behavior.
    const bool donorMajority =
        pimple.dict().getOrDefault<bool>("mcvofDonorMajority", true);
    if (donorMajority)
    {
        const scalarField& aOld = alpha1.oldTime().primitiveField();
        for (label d = 0; d < 3; ++d)
        {
            vectorField& sL = uSlopesL[d].primitiveFieldRef();
            vectorField& sG = uSlopesG[d].primitiveFieldRef();
            forAll(aOld, celli)
            {
                if (aOld[celli] < 0.5) sL[celli] = vector::zero;
                if (aOld[celli] > 0.5) sG[celli] = vector::zero;
            }
            uSlopesL[d].correctBoundaryConditions();
            uSlopesG[d].correctBoundaryConditions();
        }
    }

    const surfaceVectorField uRecL
    (
        mcvofDonorValue
        (
            mesh, alphaPhiL, U.oldTime(), uSlopesL, halfDt,
            fluxCentroid ? &tEvalLPtr() : nullptr
        )
    );
    const surfaceVectorField uRecG
    (
        mcvofDonorValue
        (
            mesh, alphaPhiG, U.oldTime(), uSlopesG, halfDt,
            fluxCentroid ? &tEvalGPtr() : nullptr
        )
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

    // spike logger (default off): report cells whose recovered
    // velocity changes by more than mcvofSpikeLogU [m/s] in one step.
    // clipResidue is alpha_new minus the alpha implied by the fluxed
    // volumes alone, so a nonzero value at an event means isoAdvector
    // corrected alpha (clip/snap) without a matching momentum
    // correction, and the mass-momentum pairing is the fault.
    {
        const scalar spikeU =
            pimple.dict().getOrDefault<scalar>("mcvofSpikeLogU", 0.0);
        if (spikeU > 0)
        {
            const volScalarField aPred
            (
                "aPred",
                alpha1.oldTime()
              - runTime.deltaT()*fvc::surfaceIntegrate(alphaPhiL)
            );
            const vectorField& uN = Ustar.primitiveField();
            const vectorField& uO = U.oldTime().primitiveField();
            const scalarField& a0c = alpha1.oldTime().primitiveField();
            const scalarField& a1c = alpha1.primitiveField();
            const scalarField& aPc = aPred.primitiveField();
            const scalarField& rhoc = rho.primitiveField();
            const vectorField& Cc = mesh.C().primitiveField();
            label nLog = 0;
            forAll(uN, celli)
            {
                if (mag(uN[celli] - uO[celli]) > spikeU && nLog < 20)
                {
                    Pout<< "mcvofSpike t=" << runTime.timeName()
                        << " cell=" << celli
                        << " x=" << Cc[celli]
                        << " a0=" << a0c[celli]
                        << " a1=" << a1c[celli]
                        << " clipResidue=" << (a1c[celli] - aPc[celli])
                        << " rho=" << rhoc[celli]
                        << " uOld=" << uO[celli]
                        << " uNew=" << uN[celli] << nl;
                    ++nLog;
                }
            }
        }
    }

    U = Ustar;
    U.correctBoundaryConditions();
}

// ************************************************************************* //
