/*---------------------------------------------------------------------------*\
    MCVOF explicit momentum advection (from interIsoConserving4Foam).

    Zhang, Popinet, Ling (2020), IJMF 130:103362, Eqs. 12-14 and 18-19,
    with isoAdvector's unsplit geometric liquid volume flux (unsplit
    consistency per Kuhn, Deskos, Sprague 2023, Sec. 3.2).

    Per step, after the alpha advection:
      liquid volume flux  alphaPhiL = advector.alphaPhi()      [m3/s]
      gas volume flux     alphaPhiG = phi - alphaPhiL          [m3/s]
      F_ml = rho1 * u_donor(rec) * alphaPhiL                  (Eq. 18)
      F_mg = rho2 * u_donor(rec) * alphaPhiG                  (Eq. 19)
      q*   = q^n - dt * div(F_ml + F_mg)                      (Eqs. 12+13)
      u*   = q* / rho^{n+1}                                   (Eq. 14)
    With the alpha staggering the alpha levels are the half-integer
    times, so alpha1.oldTime() below is alpha^{n-1/2} and the update
    reads
      q*  = rho(alpha^{n-1/2}) u^n - dt div(F)     (Basilisk pairing)
      u*  = q* / rho(alpha^{n+1/2})
    which matches conserving.h (f lagging u by dt/2, recovery through
    rho(f^{new})).
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

    surfaceVectorField uRecL
    (
        mcvofDonorValue(mesh, alphaPhiL, U.oldTime(), uSlopesL, halfDt)
    );
    surfaceVectorField uRecG
    (
        mcvofDonorValue(mesh, alphaPhiG, U.oldTime(), uSlopesG, halfDt)
    );

    // Favre-averaged face velocities (Kuhn, Deskos, Sprague 2023,
    // Eqs. 24/35): the face velocity of each phase flux is the ratio
    // of the donor-reconstructed phase momentum to the
    // donor-reconstructed phase content, evaluated untraced; the
    // staggering already centres the flux in time. Pairs the momentum
    // face velocity with the mass actually carried by the flux, which
    // the plain donor extrapolation of u does not. Where the
    // reconstructed phase content is below mcvofFavreTol the plain
    // donor value is kept.
    if (pimple.dict().getOrDefault<bool>("mcvofFavreFaces", false))
    {
        const scalar favreTol =
            pimple.dict().getOrDefault<scalar>("mcvofFavreTol", maskTol);

        const volVectorField qLfld(alpha1.oldTime()*U.oldTime());
        const volVectorField qGfld
        (
            (scalar(1) - alpha1.oldTime())*U.oldTime()
        );
        const volVectorField mLfld
        (
            alpha1.oldTime()*dimensionedVector(dimless, vector(1, 0, 0))
        );
        const volVectorField mGfld
        (
            (scalar(1) - alpha1.oldTime())
           *dimensionedVector(dimless, vector(1, 0, 0))
        );

        PtrList<volVectorField> sQL
        (
            mcvofSlopes(mesh, qLfld, limiterTheta,
                        alpha1.oldTime(), true, maskTol, false)
        );
        PtrList<volVectorField> sML
        (
            mcvofSlopes(mesh, mLfld, limiterTheta,
                        alpha1.oldTime(), true, maskTol, false)
        );
        PtrList<volVectorField> sQG
        (
            mcvofSlopes(mesh, qGfld, limiterTheta,
                        alpha1.oldTime(), false, maskTol, false)
        );
        PtrList<volVectorField> sMG
        (
            mcvofSlopes(mesh, mGfld, limiterTheta,
                        alpha1.oldTime(), false, maskTol, false)
        );

        const surfaceVectorField qLf
        (
            mcvofDonorValue(mesh, alphaPhiL, qLfld, sQL, 0.0)
        );
        const surfaceVectorField mLf
        (
            mcvofDonorValue(mesh, alphaPhiL, mLfld, sML, 0.0)
        );
        const surfaceVectorField qGf
        (
            mcvofDonorValue(mesh, alphaPhiG, qGfld, sQG, 0.0)
        );
        const surfaceVectorField mGf
        (
            mcvofDonorValue(mesh, alphaPhiG, mGfld, sMG, 0.0)
        );

        auto favreReplace = []
        (
            surfaceVectorField& uRec,
            const surfaceVectorField& qf,
            const surfaceVectorField& mf,
            const scalar tol
        )
        {
            vectorField& u = uRec.primitiveFieldRef();
            const vectorField& q = qf.primitiveField();
            const vectorField& m = mf.primitiveField();
            forAll(u, facei)
            {
                if (m[facei].x() > tol)
                {
                    u[facei] = q[facei]/m[facei].x();
                }
            }
            forAll(uRec.boundaryField(), patchi)
            {
                vectorField& up = uRec.boundaryFieldRef()[patchi];
                const vectorField& qp = qf.boundaryField()[patchi];
                const vectorField& mp = mf.boundaryField()[patchi];
                forAll(up, facei)
                {
                    if (mp[facei].x() > tol)
                    {
                        up[facei] = qp[facei]/mp[facei].x();
                    }
                }
            }
        };

        const word favrePhases =
            pimple.dict().getOrDefault<word>("mcvofFavrePhases", "both");
        if (favrePhases == "both" || favrePhases == "liquid")
        {
            favreReplace(uRecL, qLf, mLf, favreTol);
        }
        if (favrePhases == "both" || favrePhases == "gas")
        {
            favreReplace(uRecG, qGf, mGf, favreTol);
        }
    }

    // per-phase momentum update, Zhang Eqs. 12+13, carried as two
    // fields so the wisp parts can be zeroed before the density
    // division (the splitAdvector treatment; absorbs isoAdvector's
    // clip/snap residues, which alter alpha at tolerance level
    // without compensating fluxes)
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

    q1Star -= runTime.deltaT()*fvc::surfaceIntegrate(rho1*uRecL*alphaPhiL);
    q2Star -= runTime.deltaT()*fvc::surfaceIntegrate(rho2*uRecG*alphaPhiG);

    // wisp cleanup: a phase with content below tolerance in a cell
    // carries no momentum there; the momentum removed per cell and
    // step is bounded by rho*mcvofWispTol*|u|*V
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

    // velocity recovery, Zhang Eq. 14; rho already holds the new
    // mixture density (set above in alphaEqn.H)
    Ustar.primitiveFieldRef() =
        (q1Star.primitiveField() + q2Star.primitiveField())
       /rho.primitiveField();
    Ustar.correctBoundaryConditions();

    U = Ustar;
    U.correctBoundaryConditions();
}

// ************************************************************************* //
