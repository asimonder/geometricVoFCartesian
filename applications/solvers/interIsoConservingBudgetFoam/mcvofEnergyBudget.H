/*---------------------------------------------------------------------------*\
    Per-stage discrete kinetic-energy budget for the MCVOF solver
    interIsoConservingBudgetFoam.

    Stage boundaries (one begin, five/seven marks, one write per step):
      begin  : after ++runTime, before the PIMPLE loop.   E^n snapshot.
      ALPHA  : after alphaEqnSubCycle + mixture.correct.  rho^n -> rho^{n+1}
               at frozen U (bookkeeping term).
      ADV    : after mcvofAdvection (U = Ustar).          explicit MCVOF
               momentum advection at rho^{n+1}.
      VISC   : after UEqn solve.                          implicit viscous/Sp
               momentum predictor (contains the gravity, surface-tension and
               old-pressure reconstruct RHS of UEqn.H).
      P1..P3 : after each pEqn corrector.                 projection change
               per corrector (contains phig placement, pEqn.H).
      write  : after the PIMPLE loop.                     E^{n+1} snapshot,
               closing residual, physics integrals, one ASCII row.

    The stage deltas telescope, so
      dE_step = dE_alpha + dE_adv + dE_visc + dE_p1 + dE_p2 + dE_p3 + resid
    holds by construction; resid collects float rounding plus anything that
    moves U after the last recorded mark (e.g. turbulence->correct, or the
    frozenFlow early-continue which skips ADV..P3).

    Energies (Joules, global sums):
      Etot = sum 0.5*rho*|U|^2*V           rho = alpha1*rho1 + alpha2*rho2
      Ewat = sum 0.5*rho1*alpha1*|U|^2*V   water-weighted
    Physics integrals per step (Joules, evaluated on end-of-step fields):
      diss  = dt * sum 2*mu*|symm(grad U)|^2 * V     mu = mixture.mu()
      Pgrav = dt * sum rho*(g & U) * V
      Pprgh = -dt * sum U & grad(p_rgh) * V
    so the numerical residual of any stage can be separated from resolved
    physics in postprocessing.

    Activation: PIMPLE switch  mcvofEnergyBudget  (default false). Every
    hook returns immediately when the switch is false; no fields are
    allocated and no file is opened, so the inactive overhead is one bool
    test per hook.

    Output: energyBudget.dat in the case root, one ASCII row per step,
    std::ofstream append with .bak backup on restart, mirroring the
    hamiltonian.dat monitor convention. Column map in the file header.

    Assumes the v2 solver contract: static mesh, nOuterCorrectors 1,
    nAlphaSubCycles 1. nCorrectors > 3 is refused at startup.
\*---------------------------------------------------------------------------*/

#ifndef mcvofEnergyBudget_H
#define mcvofEnergyBudget_H

#include "fvCFD.H"
#include "uniformDimensionedFields.H"
#include <fstream>
#include <iomanip>

namespace Foam
{

class mcvofEnergyBudget
{
public:

    enum stageType
    {
        ALPHA = 0,   // alpha/rho update at frozen U
        ADV   = 1,   // explicit MCVOF momentum advection
        VISC  = 2,   // implicit viscous/Sp predictor solve
        P1    = 3,   // pressure corrector 1
        P2    = 4,   // pressure corrector 2
        P3    = 5,   // pressure corrector 3
        NSTAGES = 6
    };


private:

    // Configuration
    const bool active_;
    const label nCorr_;

    // Field references
    const fvMesh& mesh_;
    const Time& runTime_;
    const volVectorField& U_;
    const volScalarField& rho_;
    const volScalarField& alpha1_;
    const volScalarField& p_rgh_;
    const dimensionedScalar rho1_;
    const uniformDimensionedVectorField& g_;
    const immiscibleIncompressibleTwoPhaseGenericMixture& mixture_;

    // Per-step state
    scalar Etot0_, Ewat0_;           // step-start snapshot
    scalar EtotPrev_, EwatPrev_;     // last snapshot
    FixedList<scalar, NSTAGES> dTot_, dWat_;
    label pCorr_;
    scalar eBound_, vBound_;         // advection-block extras (col 36-37)

    // Spatial partition of the advection-block KE change (col 38-43):
    // interface band = cells within nBandLayers_ face-hops of any
    // internal face straddling alpha = 0.5 (band does not propagate
    // across processor boundaries; adequate for attribution)
    const label nBandLayers_;
    scalarField keCell0_;            // per-cell KE at step start [J]
    boolList inBand_;
    scalar dBand_, dWatBulk_, dAirBulk_;
    label nBand_;

    //- Rebuild the interface band from the current alpha
    void buildBand()
    {
        const labelUList& own = mesh_.owner();
        const labelUList& nei = mesh_.neighbour();
        const scalarField& a = alpha1_.primitiveField();

        inBand_.setSize(mesh_.nCells());
        inBand_ = false;
        for (label fi = 0; fi < mesh_.nInternalFaces(); ++fi)
        {
            if ((a[own[fi]] - 0.5)*(a[nei[fi]] - 0.5) < 0)
            {
                inBand_[own[fi]] = true;
                inBand_[nei[fi]] = true;
            }
        }
        boolList cur(inBand_);
        for (label layer = 1; layer < nBandLayers_; ++layer)
        {
            boolList next(mesh_.nCells(), false);
            for (label fi = 0; fi < mesh_.nInternalFaces(); ++fi)
            {
                const label o = own[fi];
                const label n = nei[fi];
                if (cur[o] && !inBand_[n]) { inBand_[n] = true; next[n] = true; }
                if (cur[n] && !inBand_[o]) { inBand_[o] = true; next[o] = true; }
            }
            cur = next;
        }
        nBand_ = 0;
        forAll(inBand_, c) if (inBand_[c]) ++nBand_;
        reduce(nBand_, sumOp<label>());
    }

    // Accumulated sums (Joules since simulation start of this run)
    FixedList<scalar, NSTAGES> cumTot_, cumWat_;
    scalar cumResTot_, cumResWat_;

    // Output
    const word fileName_;
    bool headerDone_;


    //- Global KE sums (all ranks participate; reduced values returned)
    void snapshot(scalar& Etot, scalar& Ewat) const
    {
        Etot = fvc::domainIntegrate(0.5*rho_*magSqr(U_)).value();
        Ewat = fvc::domainIntegrate(0.5*rho1_*alpha1_*magSqr(U_)).value();
    }

    //- Master-side append stream with restart backup, as in the
    //  hamiltonian.dat coded functionObject of the production case
    void writeRow(const std::string& row, const std::string& header)
    {
        if (!Pstream::master()) return;

        static bool firstCall = true;

        std::ifstream check(fileName_.c_str());
        const bool fileExists = check.good();
        check.close();

        if (firstCall && fileExists)
        {
            std::ifstream src(fileName_.c_str(), std::ios::binary);
            std::ofstream dst((fileName_ + ".bak").c_str(),
                              std::ios::binary);
            dst << src.rdbuf();
        }
        firstCall = false;

        std::ofstream file(fileName_.c_str(),
                           std::ios::out | std::ios::app);
        if (!fileExists)
        {
            file << header;
        }
        file << row;
        file.close();
    }


public:

    mcvofEnergyBudget
    (
        const pimpleControl& pimple,
        const fvMesh& mesh,
        const Time& runTime,
        const volVectorField& U,
        const volScalarField& rho,
        const volScalarField& alpha1,
        const volScalarField& p_rgh,
        const dimensionedScalar& rho1,
        const uniformDimensionedVectorField& g,
        const immiscibleIncompressibleTwoPhaseGenericMixture& mixture
    )
    :
        active_
        (
            pimple.dict().getOrDefault<bool>("mcvofEnergyBudget", false)
        ),
        nCorr_(pimple.dict().getOrDefault<label>("nCorrectors", 1)),
        mesh_(mesh),
        runTime_(runTime),
        U_(U),
        rho_(rho),
        alpha1_(alpha1),
        p_rgh_(p_rgh),
        rho1_(rho1),
        g_(g),
        mixture_(mixture),
        Etot0_(0), Ewat0_(0),
        EtotPrev_(0), EwatPrev_(0),
        dTot_(scalar(0)), dWat_(scalar(0)),
        pCorr_(0),
        eBound_(0), vBound_(0),
        nBandLayers_
        (
            pimple.dict().getOrDefault<label>("mcvofBandLayers", 3)
        ),
        keCell0_(),
        inBand_(),
        dBand_(0), dWatBulk_(0), dAirBulk_(0),
        nBand_(0),
        cumTot_(scalar(0)), cumWat_(scalar(0)),
        cumResTot_(0), cumResWat_(0),
        fileName_("energyBudget.dat"),
        headerDone_(false)
    {
        // Startup echo of the effective configuration (smoke-test rule:
        // every binary echoes what it is actually going to do).
        Info<< "mcvofEnergyBudget: active " << (active_ ? "yes" : "no");
        if (active_)
        {
            Info<< ", file " << fileName_
                << ", nCorrectors " << nCorr_
                << ", stages alpha/adv/visc/p1-p" << nCorr_ << "/resid"
                << ", rho1 " << rho1_.value()
                << ", |g| " << mag(g_.value());
        }
        Info<< nl << endl;

        if (active_ && nCorr_ > 3)
        {
            FatalErrorInFunction
                << "mcvofEnergyBudget supports nCorrectors <= 3, got "
                << nCorr_
                << exit(FatalError);
        }
    }


    //- Step-start snapshot; call once per time step before the PIMPLE loop
    void begin()
    {
        if (!active_) return;

        snapshot(Etot0_, Ewat0_);
        EtotPrev_ = Etot0_;
        EwatPrev_ = Ewat0_;
        forAll(dTot_, s)
        {
            dTot_[s] = 0;
            dWat_[s] = 0;
        }
        pCorr_ = 0;
        eBound_ = 0;
        vBound_ = 0;

        // per-cell KE at step start, for the block spatial partition
        keCell0_ =
            0.5*rho_.primitiveField()
           *magSqr(U_.primitiveField())
           *mesh_.V().field();
        dBand_ = dWatBulk_ = dAirBulk_ = 0;
    }


    //- Budget active? (guards the extra diagnostics in mcvofAdvection.H)
    bool on() const
    {
        return active_;
    }


    //- Advection-block extras: eB = KE change due to isoAdvector's
    //  bounding edit of alpha (recovery with bounded vs pure-flux
    //  density); vB = L1 volume of the bounding edit [m3/step]
    void setAdvExtras(const scalar eB, const scalar vB)
    {
        if (!active_) return;
        eBound_ = eB;
        vBound_ = vB;
    }


    //- Record the KE change since the previous snapshot into stage s
    void mark(const stageType s)
    {
        if (!active_) return;

        scalar Etot, Ewat;
        snapshot(Etot, Ewat);
        dTot_[s] += Etot - EtotPrev_;
        dWat_[s] += Ewat - EwatPrev_;
        EtotPrev_ = Etot;
        EwatPrev_ = Ewat;

        // At the end of the advection block (ADV), partition the
        // block's per-cell KE change (step start -> here, i.e. the
        // ALPHA + ADV stages combined) into interface band, water
        // bulk (alpha >= 0.5 outside the band) and air bulk.
        if (s == ADV)
        {
            buildBand();

            const scalarField& a = alpha1_.primitiveField();
            const scalarField keNow
            (
                0.5*rho_.primitiveField()
               *magSqr(U_.primitiveField())
               *mesh_.V().field()
            );

            scalar db = 0, dw = 0, da = 0;
            forAll(keNow, c)
            {
                const scalar dk = keNow[c] - keCell0_[c];
                if (inBand_[c])       db += dk;
                else if (a[c] >= 0.5) dw += dk;
                else                  da += dk;
            }
            reduce(db, sumOp<scalar>());
            reduce(dw, sumOp<scalar>());
            reduce(da, sumOp<scalar>());
            dBand_ = db;
            dWatBulk_ = dw;
            dAirBulk_ = da;
        }
    }


    //- Record one pressure corrector (call after each pEqn.H)
    void markPressure()
    {
        if (!active_) return;

        if (pCorr_ >= 3)
        {
            FatalErrorInFunction
                << "more than 3 pressure correctors marked"
                << exit(FatalError);
        }
        mark(static_cast<stageType>(P1 + pCorr_));
        ++pCorr_;
    }


    //- End-of-step: closing residual, physics integrals, ASCII row
    void write()
    {
        if (!active_) return;

        scalar Etot, Ewat;
        snapshot(Etot, Ewat);

        scalar sumTot = 0, sumWat = 0;
        forAll(dTot_, s)
        {
            sumTot += dTot_[s];
            sumWat += dWat_[s];
        }
        const scalar resTot = (Etot - Etot0_) - sumTot;
        const scalar resWat = (Ewat - Ewat0_) - sumWat;

        forAll(dTot_, s)
        {
            cumTot_[s] += dTot_[s];
            cumWat_[s] += dWat_[s];
        }
        cumResTot_ += resTot;
        cumResWat_ += resWat;

        // Physics integrals on end-of-step fields, converted to J/step
        const scalar dt = runTime_.deltaTValue();

        const volScalarField muMix(mixture_.mu());
        const volScalarField dissDens
        (
            2.0*muMix*magSqr(symm(fvc::grad(U_)))
        );
        const scalar diss =
            dt*fvc::domainIntegrate(dissDens).value();

        // resolved viscous dissipation split over the same band [J/step]
        scalar dissBand = 0, dissBulk = 0;
        {
            const scalarField dd
            (
                dissDens.primitiveField()*mesh_.V().field()
            );
            forAll(dd, c)
            {
                if (inBand_.size() && inBand_[c]) dissBand += dd[c];
                else dissBulk += dd[c];
            }
            reduce(dissBand, sumOp<scalar>());
            reduce(dissBulk, sumOp<scalar>());
            dissBand *= dt;
            dissBulk *= dt;
        }

        const scalar Pgrav =
            dt*fvc::domainIntegrate(rho_*(g_ & U_)).value();

        const scalar Pprgh =
           -dt*fvc::domainIntegrate(U_ & fvc::grad(p_rgh_)).value();

        // One ASCII row
        std::ostringstream os;
        os << std::setprecision(16) << std::scientific;
        os << runTime_.timeOutputValue() << " " << dt;
        forAll(dTot_, s) os << " " << dTot_[s];
        os << " " << resTot;
        forAll(dWat_, s) os << " " << dWat_[s];
        os << " " << resWat;
        os << " " << Etot << " " << Ewat;
        os << " " << diss << " " << Pgrav << " " << Pprgh;
        forAll(cumTot_, s) os << " " << cumTot_[s];
        os << " " << cumResTot_;
        forAll(cumWat_, s) os << " " << cumWat_[s];
        os << " " << cumResWat_;
        os << " " << eBound_ << " " << vBound_;
        os << " " << dBand_ << " " << dWatBulk_ << " " << dAirBulk_;
        os << " " << dissBand << " " << dissBulk << " " << nBand_;
        os << "\n";

        static const std::string header =
            "# mcvofEnergyBudget: KE budget per time step [J]\n"
            "# 1:t 2:dt"
            " 3:dEtot_alpha 4:dEtot_adv 5:dEtot_visc"
            " 6:dEtot_p1 7:dEtot_p2 8:dEtot_p3 9:dEtot_resid"
            " 10:dEwat_alpha 11:dEwat_adv 12:dEwat_visc"
            " 13:dEwat_p1 14:dEwat_p2 15:dEwat_p3 16:dEwat_resid"
            " 17:Etot 18:Ewat"
            " 19:diss 20:Pgrav 21:Pprgh"
            " 22-28:cumTot(alpha,adv,visc,p1,p2,p3,resid)"
            " 29-35:cumWat(alpha,adv,visc,p1,p2,p3,resid)"
            " 36:eBound_J 37:vBound_m3"
            " 38:dEblk_band 39:dEblk_watBulk 40:dEblk_airBulk"
            " 41:dissBand 42:dissBulk 43:nBandCells\n";

        writeRow(os.str(), header);

        // Terse per-step echo for the smoke run
        Info<< "energyBudget: dE " << (Etot - Etot0_)
            << " a " << dTot_[ALPHA]
            << " adv " << dTot_[ADV]
            << " visc " << dTot_[VISC]
            << " p " << dTot_[P1] + dTot_[P2] + dTot_[P3]
            << " res " << resTot << endl;
    }
};

} // End namespace Foam

#endif

// ************************************************************************* //
