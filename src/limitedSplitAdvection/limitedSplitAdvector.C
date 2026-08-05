/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2026 Asim Onder
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "limitedSplitAdvector.H"
#include "addToRunTimeSelectionTable.H"
#include "interfaceBand.H"
#include "coupledFvPatchField.H"
#include "cellSet.H"
#include "emptyFvPatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(limitedSplitAdvector, 0);
    addToRunTimeSelectionTable
    (
        geometricAdvector,
        limitedSplitAdvector,
        components
    );
}

// * * * * * * * * * * * * * * * * Constructor  * * * * * * * * * * * * * * //

Foam::limitedSplitAdvector::limitedSplitAdvector
(
    volScalarField& alpha1,
    const surfaceScalarField& phi,
    const volVectorField& U
)
:
    alpha1_(alpha1),
    phi_(phi),
    U_(U),
    mesh_(alpha1.mesh()),
    dict_(mesh_.solverDict(alpha1.name())),
    surfPtr_(reconstructionSchemes::New(alpha1_, phi_, U_, dict_)),
    alphaPhi_
    (
        IOobject
        (
            "wyAlphaPhi",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(phi.dimensions(), Zero)
    ),
    qPhi_
    (
        IOobject
        (
            "wyQPhi",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector
        (
            phi.dimensions()*dimDensity*dimVelocity, Zero
        )
    ),
    sweepDirs_(),
    nSweeps_(0),
    inZone_(),
    inMargin_(),
    wholeDomain_(true),
    dynamicBand_(false),
    nCellsAtZoneBuild_(0),
    cH_(mesh_.nCells(), 0),
    clampTol_(dict_.getOrDefault<scalar>("wyClampTol", 1e-12)),
    fatalTol_(dict_.getOrDefault<scalar>("wyFatalTol", 1e-9)),
    slopeLimiter_
    (
        dict_.getOrDefault<word>("slopeLimiter", "none") == "minmod" ? 1
      : dict_.getOrDefault<word>("slopeLimiter", "none") == "vanLeer" ? 2
      : 0
    ),
    limitBandOnly_
    (
        dict_.getOrDefault<word>("slopeLimiterScope", "band") == "band"
    ),
    nClamped_(0),
    maxExcursion_(0)
{
    {
        const word lim(dict_.getOrDefault<word>("slopeLimiter", "none"));
        const word scope(dict_.getOrDefault<word>("slopeLimiterScope", "band"));
        if (scope != "band" && scope != "domain")
        {
            FatalErrorInFunction
                << "slopeLimiterScope must be band or domain, got "
                << scope << exit(FatalError);
        }
        if (lim != "none" && lim != "minmod" && lim != "vanLeer")
        {
            FatalErrorInFunction
                << "slopeLimiter must be none, minmod or vanLeer, got "
                << lim << exit(FatalError);
        }
        Info<< "limitedSplitAdvector: majority-phase slope limiter "
            << lim << ", scope " << scope << endl;
    }
    // --- Scope guards (stage 1) ---------------------------------------

    // exactly two active directions (2D with one empty direction)
    const Vector<label>& geoD = mesh_.geometricD();
    nSweeps_ = 0;
    for (direction c = 0; c < 3; ++c)
    {
        sweepDirs_[c] = -1;
        if (geoD[c] == 1)
        {
            sweepDirs_[nSweeps_++] = c;
        }
    }
    if (nSweeps_ < 2)
    {
        FatalErrorInFunction
            << "limitedSplitAdvector needs 2 or 3 active directions; geometricD = "
            << geoD << exit(FatalError);
    }

    buildZone();

    checkZone();

    if (dict_.found("wyMomentumBand") || dict_.found("wyMomentumBandLayers"))
    {
        FatalErrorInFunction
            << "wyMomentumBand has been removed: momentum is advected "
            << "on all cells by the per-phase Lax-Wendroff sweep "
            << "fluxes. Remove the entries from fvSolution."
            << exit(FatalError);
    }

    word dirsStr =
        Foam::name(sweepDirs_[0]) + " " + Foam::name(sweepDirs_[1]);
    if (nSweeps_ == 3)
    {
        dirsStr += " " + Foam::name(sweepDirs_[2]);
    }
    Info<< "limitedSplitAdvector: " << nSweeps_ << " sweeps ("
        << dirsStr
        << ") alternating, reconstruction "
        << dict_.getOrDefault<word>("reconstructionScheme", "plicRDF")
        << ", clampTol " << clampTol_
        << ", fatalTol " << fatalTol_
        << ", cells " << returnReduce(mesh_.nCells(), sumOp<label>());
    if (wholeDomain_)
    {
        Info<< ", zone domain";
    }
    else
    {
        label nZ = 0;
        forAll(inZone_, c) if (inZone_[c]) ++nZ;
        Info<< ", zone " << returnReduce(nZ, sumOp<label>())
            << " cells";
    }
    Info<< nl << endl;
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::limitedSplitAdvector::checkZone()
{
    // uniform cell volume over the zone
    const scalarField& V = mesh_.V().field();
    scalar V0 = -1;
    forAll(V, celli)
    {
        if (inZone_[celli]) { V0 = V[celli]; break; }
    }
    reduce(V0, maxOp<scalar>());
    if (V0 < 0)
    {
        FatalErrorInFunction
            << "limitedSplitAdvector zone is empty" << exit(FatalError);
    }
    scalar dVmax = 0;
    forAll(V, celli)
    {
        if (inZone_[celli])
        {
            dVmax = max(dVmax, mag(V[celli] - V0));
        }
    }
    reduce(dVmax, maxOp<scalar>());
    if (dVmax > 1e-10*V0)
    {
        FatalErrorInFunction
            << "limitedSplitAdvector requires a uniform mesh over its zone; "
            << "max relative cell volume deviation = " << dVmax/V0
            << exit(FatalError);
    }

    // axis-aligned faces, domain-wide: the momentum sweeps visit
    // every face, not only the lattice
    const labelUList& ownc = mesh_.owner();
    const surfaceVectorField& Sf = mesh_.Sf();
    const surfaceScalarField& magSf = mesh_.magSf();
    scalar worstAlign = 1;
    forAll(ownc, facei)
    {
        const vector nHat = Sf[facei]/magSf[facei];
        worstAlign = min(worstAlign, cmptMax(cmptMag(nHat)));
    }
    reduce(worstAlign, minOp<scalar>());
    if (worstAlign < 1 - 1e-10)
    {
        FatalErrorInFunction
            << "limitedSplitAdvector requires axis-aligned faces (momentum "
            << "sweeps cover the whole mesh); worst alignment = "
            << worstAlign << exit(FatalError);
    }
}


void Foam::limitedSplitAdvector::buildZone()
{
    inZone_.setSize(mesh_.nCells());
    inZone_ = true;
    wholeDomain_ = true;

    IOobject dictIO
    (
        "ijkDict",
        mesh_.time().system(),
        mesh_,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    );

    word setName("domain");
    dynamicBand_ = false;
    if (dictIO.typeHeaderOk<IOdictionary>(true))
    {
        IOdictionary ijkDict(dictIO);
        setName = ijkDict.getOrDefault<word>("cellSet", "domain");
        dynamicBand_ = ijkDict.getOrDefault<bool>("dynamicBand", false);
    }

    if (dynamicBand_)
    {
        // finest-level band: every cell at the minimum cell volume;
        // the margin logic below then marks refinement-level jumps
        wholeDomain_ = false;
        const scalarField& V = mesh_.V().field();
        const scalar Vmin = returnReduce(min(V), minOp<scalar>());
        inZone_ = false;
        forAll(V, celli)
        {
            if (V[celli] < (1 + 1e-6)*Vmin) inZone_[celli] = true;
        }
    }
    else if (setName != "domain")
    {
        wholeDomain_ = false;
        inZone_ = false;
        cellSet zoneCells(mesh_, setName);
        for (const label celli : zoneCells)
        {
            inZone_[celli] = true;
        }
    }

    // confinement margin: zone cells within two face-hops of any
    // non-zone cell (mixed alpha here is fatal at every step)
    inMargin_.setSize(mesh_.nCells());
    inMargin_ = false;
    if (!wholeDomain_)
    {
        const labelUList& own = mesh_.owner();
        const labelUList& nei = mesh_.neighbour();
        boolList edge(mesh_.nCells(), false);
        forAll(own, facei)
        {
            const bool o = inZone_[own[facei]];
            const bool n = inZone_[nei[facei]];
            if (o != n)
            {
                if (o) edge[own[facei]] = true;
                if (n) edge[nei[facei]] = true;
            }
        }
        inMargin_ = edge;
        boolList cur(edge);
        for (label layer = 1; layer < 2; ++layer)
        {
            boolList next(mesh_.nCells(), false);
            forAll(own, facei)
            {
                const label o = own[facei];
                const label n = nei[facei];
                if (cur[o] && inZone_[n] && !inMargin_[n])
                {
                    inMargin_[n] = true; next[n] = true;
                }
                if (cur[n] && inZone_[o] && !inMargin_[o])
                {
                    inMargin_[o] = true; next[o] = true;
                }
            }
            cur = next;
        }
    }

    nCellsAtZoneBuild_ = mesh_.nCells();
}


Foam::scalar Foam::limitedSplitAdvector::slabFraction
(
    const scalar mx,
    const scalar mt,
    const scalar a,
    const scalar xa,
    const scalar xb
)
{
    // Unit cell, liquid region mx*x + mt*t <= ai, mx,mt >= 0,
    // mx + mt = 1, total liquid fraction a in (0,1). Returns the
    // liquid fraction of the slab [xa,xb] x [0,1].

    const scalar w = xb - xa;
    if (w <= SMALL) return a;

    const scalar eps = 1e-12;

    // interface parallel to the sweep direction: uniform column height
    if (mx < eps)
    {
        return a;                       // h(x) = a everywhere
    }

    // interface normal to the sweep direction: sharp front at x = a
    if (mt < eps)
    {
        // liquid occupies x <= ai with ai = a
        const scalar xl = max(xa, min(xb, a));
        return (xl - xa)/w;
    }

    // generic case: invert a -> ai
    const scalar m1 = min(mx, mt);
    const scalar m2 = max(mx, mt);
    const scalar aSym = min(a, 1 - a);
    scalar aiSym;
    if (aSym <= 0)
    {
        aiSym = 0;
    }
    else if (aSym < 0.5*m1/m2)
    {
        aiSym = sqrt(2*aSym*m1*m2);
    }
    else
    {
        aiSym = aSym*m2 + 0.5*m1;
    }
    const scalar ai = (a <= 0.5) ? aiSym : (1 - aiSym);

    // column height h(x) = clamp((ai - mx*x)/mt, 0, 1), decreasing
    // breakpoints: h = 1 at x1, h = 0 at x0, x1 <= x0
    const scalar x1 = (ai - mt)/mx;
    const scalar x0 = ai/mx;

    // integrate over [xa, xb]
    const scalar lo1 = min(max(xa, scalar(-GREAT)), xb);
    scalar area = 0;

    // full part: [xa, min(xb, x1)] with h = 1
    const scalar bFull = min(xb, x1);
    if (bFull > lo1)
    {
        area += (bFull - xa);
    }

    // linear part: [max(xa, x1), min(xb, x0)]
    const scalar aLin = max(xa, x1);
    const scalar bLin = min(xb, x0);
    if (bLin > aLin)
    {
        const scalar xm = 0.5*(aLin + bLin);
        area += (bLin - aLin)*(ai - mx*xm)/mt;
    }

    return min(max(area/w, 0.0), 1.0);
}


Foam::scalar Foam::limitedSplitAdvector::boxFraction3
(
    const scalar m1in,
    const scalar m2in,
    const scalar m3in,
    const scalar aiIn
)
{
    // sort so m1 <= m2 <= m3; L1-normalized non-negative components
    scalar m1 = mag(m1in), m2 = mag(m2in), m3 = mag(m3in);
    const scalar mSum = m1 + m2 + m3;
    if (mSum < SMALL)
    {
        return min(max(aiIn, 0.0), 1.0);
    }
    m1 /= mSum; m2 /= mSum; m3 /= mSum;
    if (m1 > m2) { const scalar t = m1; m1 = m2; m2 = t; }
    if (m2 > m3) { const scalar t = m2; m2 = m3; m3 = t; }
    if (m1 > m2) { const scalar t = m1; m1 = m2; m2 = t; }

    // symmetric wrapper: evaluate for ai <= 1/2
    const scalar ai = min(max(aiIn, 0.0), 1.0);
    const scalar a5 = min(ai, 1 - ai);
    const scalar eps = 1e-12;

    scalar F;
    if (a5 <= 0)
    {
        F = 0;
    }
    else if (m3 < eps)
    {
        F = a5;                          // degenerate: no normal
    }
    else if (m2 < eps)
    {
        F = min(a5/m3, 0.5);     // 1D front (m3 = 1)
    }
    else if (m1 < eps)
    {
        // 2D wedge in (y,z): line m2 y + m3 z = a5
        if (a5 < m2)
        {
            F = a5*a5/(2*m2*m3);
        }
        else
        {
            F = (a5 - 0.5*m2)/m3;
        }
    }
    else
    {
        // full 3D Scardovelli-Zaleski branches, a5 <= 1/2
        if (a5 < m1)
        {
            F = a5*a5*a5/(6*m1*m2*m3);
        }
        else if (a5 < m2)
        {
            F = (3*a5*(a5 - m1) + m1*m1)/(6*m2*m3);
        }
        else if (a5 < min(m3, m1 + m2))
        {
            F =
            (
                a5*a5*a5
              - pow3(a5 - m1)
              - pow3(a5 - m2)
            )/(6*m1*m2*m3);
        }
        else if (m3 < m1 + m2)
        {
            // a5 >= m3: subtract the third corner too
            F =
            (
                a5*a5*a5
              - pow3(a5 - m1)
              - pow3(a5 - m2)
              - pow3(a5 - m3)
            )/(6*m1*m2*m3);
        }
        else
        {
            // m3 >= m1 + m2 and a5 >= m1 + m2: linear regime
            F = (2*a5 - m1 - m2)/(2*m3);
        }
    }

    return (ai <= 0.5) ? F : (1 - F);
}


Foam::scalar Foam::limitedSplitAdvector::boxIntercept3
(
    const scalar m1,
    const scalar m2,
    const scalar m3,
    const scalar a
)
{
    if (a <= 0) return 0;
    if (a >= 1) return 1;

    // pure bisection on the monotone forward function: no derivative,
    // no floating-point trap surface, 1e-15 in ~50 halvings
    scalar lo = 0, hi = 1;
    for (label it = 0; it < 52; ++it)
    {
        const scalar mid = 0.5*(lo + hi);
        if (boxFraction3(m1, m2, m3, mid) < a) lo = mid;
        else hi = mid;
    }
    return 0.5*(lo + hi);
}


Foam::scalar Foam::limitedSplitAdvector::donorFraction
(
    const label d,
    const scalar sOut,
    const scalar w,
    const scalar aD,
    const vector& nD
) const
{
    if (aD <= 0) return 0;
    if (aD >= 1) return 1;

    // interface normal pointing from liquid into gas
    const vector m = -nD;

    // slab interval in the donor's local d coordinate
    scalar xa = (sOut > 0) ? (1 - w) : 0;
    scalar xb = (sOut > 0) ? 1.0 : w;

    if (nSweeps_ == 2)
    {
        const label t =
            (sweepDirs_[0] == d) ? sweepDirs_[1] : sweepDirs_[0];
        scalar mx = m[d];
        scalar mt = m[t];
        const scalar mSum = mag(mx) + mag(mt);
        if (mSum < SMALL) return aD;
        mx /= mSum;
        mt /= mSum;
        if (mx < 0)
        {
            mx = -mx;
            const scalar t1 = 1 - xb;
            const scalar t2 = 1 - xa;
            xa = t1;
            xb = t2;
        }
        if (mt < 0) mt = -mt;   // transverse mirror: no slab effect
        return slabFraction(mx, mt, aD, xa, xb);
    }

    // --- 3D path: plane-cut cuboid sub-slab -----------------------
    const label t1i = (d + 1) % 3;
    const label t2i = (d + 2) % 3;
    scalar mx = m[d];
    scalar m1 = mag(m[t1i]);            // transverse mirrors free
    scalar m2 = mag(m[t2i]);
    const scalar mSum = mag(mx) + m1 + m2;
    if (mSum < SMALL) return aD;
    mx /= mSum; m1 /= mSum; m2 /= mSum;
    if (mx < 0)
    {
        mx = -mx;
        const scalar u1 = 1 - xb;
        const scalar u2 = 1 - xa;
        xa = u1;
        xb = u2;
    }

    // intercept from the donor's total fraction, full cell
    const scalar ai = boxIntercept3(mx, m1, m2, aD);

    // sub-slab [xa,xb]: substitute x = xa + wS*x' and renormalize
    const scalar wS = xb - xa;
    if (wS <= SMALL) return aD;
    const scalar sN = mx*wS + m1 + m2;
    const scalar aiS = (ai - mx*xa)/sN;
    return boxFraction3(mx*wS/sN, m1/sN, m2/sN, aiS);
}


template<class Type>
void Foam::limitedSplitAdvector::gatherNeighbours
(
    const label d,
    const GeometricField<Type, fvPatchField, volMesh>& f,
    Field<Type>& fMinus,
    Field<Type>& fPlus
) const
{
    const labelUList& own = mesh_.owner();
    const labelUList& nei = mesh_.neighbour();
    const surfaceVectorField& Sf = mesh_.Sf();
    const surfaceScalarField& magSf = mesh_.magSf();

    fMinus = f.primitiveField();      // physical-boundary fallback:
    fPlus = f.primitiveField();       // copy own value (zero gradient)

    forAll(own, facei)
    {
        const scalar dirComp = Sf[facei][d]/magSf[facei];
        if (mag(dirComp) < 0.99) continue;

        if (dirComp > 0)
        {
            fPlus[own[facei]] = f[nei[facei]];
            fMinus[nei[facei]] = f[own[facei]];
        }
        else
        {
            fMinus[own[facei]] = f[nei[facei]];
            fPlus[nei[facei]] = f[own[facei]];
        }
    }

    forAll(mesh_.boundary(), patchi)
    {
        const fvPatch& patch = mesh_.boundary()[patchi];
        if (!patch.coupled()) continue;

        const labelUList& fc = patch.faceCells();
        const vectorField& pSf = Sf.boundaryField()[patchi];
        const scalarField& pMagSf = magSf.boundaryField()[patchi];
        const Field<Type> pNei
        (
            refCast<const coupledFvPatchField<Type>>
            (
                f.boundaryField()[patchi]
            ).patchNeighbourField()
        );

        forAll(fc, fi)
        {
            const scalar dirComp = pSf[fi][d]/pMagSf[fi];
            if (mag(dirComp) < 0.99) continue;

            if (dirComp > 0) fPlus[fc[fi]] = pNei[fi];
            else             fMinus[fc[fi]] = pNei[fi];
        }
    }
}


Foam::vector Foam::limitedSplitAdvector::concGradient
(
    const scalar cM, const scalar cC, const scalar cP,
    const vector& tM, const vector& tC, const vector& tP,
    const bool inverse,
    const label limiter
)
{
    // Majority-phase concentration gradient (vof.h,
    // vof_concentration_gradient_x): centered where donor and both
    // neighbours carry the phase (content >= 1/2), one-sided away
    // from the interface where one neighbour drops out, zero where
    // the donor itself does not carry the phase. Gradient per unit
    // cell width (the caller works in cell-fraction coordinates).
    const scalar eps = 1e-10;
    const scalar pM = inverse ? (1 - cM) : cM;
    const scalar pC = inverse ? (1 - cC) : cC;
    const scalar pP = inverse ? (1 - cP) : cP;

    if (pC < 0.5) return vector::zero;

    const vector uM = (pM > eps) ? tM/pM : vector::zero;
    const vector uC = tC/max(pC, eps);
    const vector uP = (pP > eps) ? tP/pP : vector::zero;

    if (pP >= 0.5)
    {
        if (pM >= 0.5)
        {
            if (limiter == 0) return 0.5*(uP - uM);   // unlimited centred

            const vector dM(uC - uM);
            const vector dP(uP - uC);
            vector g(vector::zero);
            for (direction c = 0; c < 3; ++c)
            {
                const scalar a = dM[c];
                const scalar b = dP[c];
                if (a*b <= 0) continue;
                if (limiter == 1)
                {
                    // minmod
                    g[c] = (mag(a) < mag(b)) ? a : b;
                }
                else
                {
                    // van Leer harmonic
                    g[c] = 2*a*b/(a + b);
                }
            }
            return g;
        }
        return (uP - uC);                       // one-sided +
    }
    if (pM >= 0.5)
    {
        return (uC - uM);                       // one-sided -
    }
    return vector::zero;
}


Foam::scalar Foam::limitedSplitAdvector::sweep
(
    const label d,
    const scalar dt,
    volVectorField* q1,
    volVectorField* q2
)
{
    surfPtr_->reconstruct();

    // The processor/cyclic halo of the reconstruction normal must be
    // fresh: both sides of a coupled face compute the donor flux
    // independently, and conservation requires them to compute it
    // from bitwise-identical inputs. A stale normal halo broke mass
    // conservation on decomposed meshes (measured 0.46 percent per
    // step at 4 ranks before this call).
    surfPtr_->normal().correctBoundaryConditions();

    const labelUList& own = mesh_.owner();
    const labelUList& nei = mesh_.neighbour();
    const surfaceVectorField& Sf = mesh_.Sf();
    const surfaceScalarField& magSf = mesh_.magSf();
    const scalarField& V = mesh_.V().field();
    const volVectorField& normal = surfPtr_->normal();
    const scalarField& a = alpha1_.primitiveField();
    const scalar aTol = 1e-8;   // face-alpha band: only interface faces feed the WY 0.5 check

    scalarField divF(mesh_.nCells(), 0);      // [m3/s]
    scalarField divPhiD(mesh_.nCells(), 0);   // [m3/s]
    scalar maxCo = 0;

    // stage 2: momentum sweep data
    const bool withQ = (q1 != nullptr);

    vectorField divFq1, divFq2;
    autoPtr<volVectorField> g1Fld, g2Fld;
    if (withQ)
    {
        divFq1.setSize(mesh_.nCells(), vector::zero);
        divFq2.setSize(mesh_.nCells(), vector::zero);

        scalarField alM(mesh_.nCells()), alP(mesh_.nCells());
        vectorField q1M(mesh_.nCells()), q1P(mesh_.nCells());
        vectorField q2M(mesh_.nCells()), q2P(mesh_.nCells());
        gatherNeighbours(d, alpha1_, alM, alP);
        gatherNeighbours(d, *q1, q1M, q1P);
        gatherNeighbours(d, *q2, q2M, q2P);

        // per-sweep concentration gradients as registered fields so
        // their remote (cyclic/processor) donor values are available
        // through the coupled-patch exchange; this keeps the scheme
        // fully second order across periodic seams
        g1Fld.reset
        (
            new volVectorField
            (
                IOobject
                (
                    "wyG1", mesh_.time().timeName(), mesh_,
                    IOobject::NO_READ, IOobject::NO_WRITE,
                    IOobject::NO_REGISTER
                ),
                mesh_,
                dimensionedVector(q1->dimensions(), Zero)
            )
        );
        g2Fld.reset
        (
            new volVectorField
            (
                IOobject
                (
                    "wyG2", mesh_.time().timeName(), mesh_,
                    IOobject::NO_READ, IOobject::NO_WRITE,
                    IOobject::NO_REGISTER
                ),
                mesh_,
                dimensionedVector(q2->dimensions(), Zero)
            )
        );

        const scalarField& aIn = alpha1_.primitiveField();
        const vectorField& q1In = q1->primitiveField();
        const vectorField& q2In = q2->primitiveField();
        vectorField& g1 = g1Fld->primitiveFieldRef();
        vectorField& g2 = g2Fld->primitiveFieldRef();
        // the limiter acts inside the interface band only, so the
        // bulk turbulence keeps the unlimited centred slope
        const scalarField* bandPtr = nullptr;
        if (slopeLimiter_ != 0 && limitBandOnly_)
        {
            bandPtr =
                &interfaceBand::New(mesh_).band().primitiveField();
        }

        forAll(g1, celli)
        {
            const label lim =
                (bandPtr && (*bandPtr)[celli] < 0.5) ? 0 : slopeLimiter_;
            g1[celli] = concGradient
            (
                alM[celli], aIn[celli], alP[celli],
                q1M[celli], q1In[celli], q1P[celli],
                false, lim
            );
            g2[celli] = concGradient
            (
                alM[celli], aIn[celli], alP[celli],
                q2M[celli], q2In[celli], q2P[celli],
                true, lim
            );
        }
        g1Fld->correctBoundaryConditions();
        g2Fld->correctBoundaryConditions();
    }

    // per-phase face value at the donor's outflow face:
    // concentration + (1-Co) half-cell offset along the gradient
    auto faceValue = []
    (
        const scalar sOut, const scalar w,
        const scalar pD, const vector& qD, const vector& gD
    ) -> vector
    {
        if (pD < 1e-10) return vector::zero;
        return qD/pD + sOut*0.5*(1 - w)*gD;
    };

    // --- internal faces ------------------------------------------------
    // Alpha transport lives on the lattice (zone faces). Momentum is
    // transported on EVERY axis-aligned face: outside the lattice the
    // phases are pure, donorFraction is trivially 0 or 1, and the
    // per-phase flux below reduces to the Lax-Wendroff flux of the
    // bulk momentum. The Weymouth-Yue compression bookkeeping
    // (divPhiD) is accumulated into lattice cells only.
    forAll(own, facei)
    {
        const bool zoneFace =
            inZone_[own[facei]] || inZone_[nei[facei]];
        if (!zoneFace && !withQ) continue;

        const scalar dirComp = Sf[facei][d]/magSf[facei];
        if (mag(dirComp) < 0.99) continue;

        const scalar phiF = phi_[facei];
        const scalar dirSign = sign(dirComp);

        if (inZone_[own[facei]]) divPhiD[own[facei]] += phiF;
        if (inZone_[nei[facei]]) divPhiD[nei[facei]] -= phiF;

        if (mag(phiF) < VSMALL) continue;

        const label dc = (phiF >= 0) ? own[facei] : nei[facei];
        const scalar sOut = (dc == own[facei]) ? dirSign : -dirSign;
        const scalar w = mag(phiF)*dt/V[dc];
        const scalar af = 0.5*(a[own[facei]] + a[nei[facei]]);
        if (af > aTol && af < 1 - aTol) maxCo = max(maxCo, w);

        const scalar frac = donorFraction(d, sOut, w, a[dc], normal[dc]);
        const scalar F = frac*phiF;

        if (zoneFace)
        {
            divF[own[facei]] += F;
            divF[nei[facei]] -= F;
            alphaPhi_.primitiveFieldRef()[facei] += F;
        }

        if (withQ)
        {
            const vector ff1 = faceValue
            (
                sOut, w, a[dc], (*q1)[dc], (*g1Fld)[dc]
            );
            const vector ff2 = faceValue
            (
                sOut, w, 1 - a[dc], (*q2)[dc], (*g2Fld)[dc]
            );
            const vector Fq1 = ff1*F;
            const vector Fq2 = ff2*(phiF - F);

            divFq1[own[facei]] += Fq1;
            divFq1[nei[facei]] -= Fq1;
            divFq2[own[facei]] += Fq2;
            divFq2[nei[facei]] -= Fq2;
            qPhi_.primitiveFieldRef()[facei] += Fq1 + Fq2;
        }
    }

    // --- boundary faces ------------------------------------------------
    forAll(mesh_.boundary(), patchi)
    {
        const fvPatch& patch = mesh_.boundary()[patchi];
        if (isA<emptyFvPatch>(patch)) continue;

        const labelUList& fc = patch.faceCells();
        const vectorField& pSf = Sf.boundaryField()[patchi];
        const scalarField& pMagSf = magSf.boundaryField()[patchi];
        const scalarField& pPhi = phi_.boundaryField()[patchi];
        scalarField& pAlphaPhi = alphaPhi_.boundaryFieldRef()[patchi];

        const bool coupled = patch.coupled();

        tmp<scalarField> tANei;
        tmp<vectorField> tNNei;
        if (coupled)
        {
            tANei = refCast<const coupledFvPatchField<scalar>>
            (
                alpha1_.boundaryField()[patchi]
            ).patchNeighbourField();
            tNNei = refCast<const coupledFvPatchField<vector>>
            (
                normal.boundaryField()[patchi]
            ).patchNeighbourField();
        }

        forAll(fc, fi)
        {
            const bool zoneCell = inZone_[fc[fi]];
            if (!zoneCell && !withQ) continue;

            const scalar dirComp = pSf[fi][d]/pMagSf[fi];
            if (mag(dirComp) < 0.99) continue;

            const scalar phiF = pPhi[fi];
            const scalar dirSign = sign(dirComp);

            if (zoneCell) divPhiD[fc[fi]] += phiF;

            if (mag(phiF) < VSMALL) continue;

            scalar F = 0;
            if (phiF >= 0)
            {
                // outflow: donor is the local cell
                const scalar w = phiF*dt/V[fc[fi]];
                const scalar af =
                    0.5*(a[fc[fi]] + (coupled ? tANei()[fi] : a[fc[fi]]));
                if (af > aTol && af < 1 - aTol) maxCo = max(maxCo, w);
                F = phiF*donorFraction
                (
                    d, dirSign, w, a[fc[fi]], normal[fc[fi]]
                );
            }
            else if (coupled)
            {
                // inflow through a coupled patch: donor is remote;
                // remote volume equals local (uniform mesh, verified)
                const scalar w = mag(phiF)*dt/V[fc[fi]];
                const scalar af = 0.5*(a[fc[fi]] + tANei()[fi]);
                if (af > aTol && af < 1 - aTol) maxCo = max(maxCo, w);
                F = phiF*donorFraction
                (
                    d, -dirSign, w, tANei()[fi], tNNei()[fi]
                );
            }
            else
            {
                // inflow through a physical boundary: boundary value
                F = phiF*alpha1_.boundaryField()[patchi][fi];
            }

            if (zoneCell)
            {
                divF[fc[fi]] += F;
                pAlphaPhi[fi] += F;
            }

            if (withQ)
            {
                // outflow: local donor. coupled inflow: remote donor
                // with the remote gradient through the coupled-patch
                // exchange, second order across periodic seams.
                vector Fq1(Zero), Fq2(Zero);
                if (phiF >= 0)
                {
                    const label dc = fc[fi];
                    const scalar w2 = phiF*dt/V[dc];
                    Fq1 = F*faceValue
                    (
                        dirSign, w2, a[dc], (*q1)[dc], (*g1Fld)[dc]
                    );
                    Fq2 = (phiF - F)*faceValue
                    (
                        dirSign, w2, 1 - a[dc], (*q2)[dc], (*g2Fld)[dc]
                    );
                }
                else if (coupled)
                {
                    const scalar aN = tANei()[fi];
                    const scalar w2 = mag(phiF)*dt/V[fc[fi]];
                    const vector q1N =
                        refCast<const coupledFvPatchField<vector>>
                        (
                            q1->boundaryField()[patchi]
                        ).patchNeighbourField()()[fi];
                    const vector q2N =
                        refCast<const coupledFvPatchField<vector>>
                        (
                            q2->boundaryField()[patchi]
                        ).patchNeighbourField()()[fi];
                    const vector g1N =
                        refCast<const coupledFvPatchField<vector>>
                        (
                            g1Fld->boundaryField()[patchi]
                        ).patchNeighbourField()()[fi];
                    const vector g2N =
                        refCast<const coupledFvPatchField<vector>>
                        (
                            g2Fld->boundaryField()[patchi]
                        ).patchNeighbourField()()[fi];

                    Fq1 = F*faceValue
                    (
                        -dirSign, w2, aN, q1N, g1N
                    );
                    Fq2 = (phiF - F)*faceValue
                    (
                        -dirSign, w2, 1 - aN, q2N, g2N
                    );
                }

                divFq1[fc[fi]] += Fq1;
                divFq2[fc[fi]] += Fq2;
                qPhi_.boundaryFieldRef()[patchi][fi] += Fq1 + Fq2;
            }
        }
    }

    // --- cell update: flux + Weymouth-Yue compression -------------------
    scalarField& aRef = alpha1_.primitiveFieldRef();
    forAll(aRef, celli)
    {
        if (!inZone_[celli]) continue;
        aRef[celli] +=
            (dt/V[celli])*(-divF[celli] + cH_[celli]*divPhiD[celli]);
    }
    alpha1_.correctBoundaryConditions();

    if (withQ)
    {
        vectorField& q1Ref = q1->primitiveFieldRef();
        vectorField& q2Ref = q2->primitiveFieldRef();
        forAll(q1Ref, celli)
        {
            const scalar f = dt/V[celli];
            q1Ref[celli] +=
                f*(-divFq1[celli] + tc1_[celli]*divPhiD[celli]);
            q2Ref[celli] +=
                f*(-divFq2[celli] + tc2_[celli]*divPhiD[celli]);
        }
        q1->correctBoundaryConditions();
        q2->correctBoundaryConditions();
    }

    reduce(maxCo, maxOp<scalar>());
    return maxCo;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::limitedSplitAdvector::advect()
{
    advanceStep(nullptr, nullptr);
}


void Foam::limitedSplitAdvector::advect(volVectorField& q1, volVectorField& q2)
{
    advanceStep(&q1, &q2);
}


void Foam::limitedSplitAdvector::advanceStep
(
    volVectorField* q1p,
    volVectorField* q2p
)
{
    const scalar dt = mesh_.time().deltaTValue();

    // AMR: rebuild the zone and re-verify its guards when the mesh
    // has changed; static meshes never trigger this
    const bool meshChanged = returnReduce
    (
        mesh_.nCells() != nCellsAtZoneBuild_, orOp<bool>()
    );
    if (meshChanged)
    {
        Info<< "limitedSplitAdvector: mesh change ("
            << nCellsAtZoneBuild_ << " -> " << mesh_.nCells()
            << " cells); rebuilding zone" << nl;
        buildZone();
        checkZone();
        cH_.setSize(mesh_.nCells(), 0);
    }

    // frozen Heaviside marker of the step (Weymouth-Yue Eq. 13)
    const scalarField& a0 = alpha1_.primitiveField();
    forAll(a0, celli)
    {
        cH_[celli] = (a0[celli] > 0.5) ? 1 : 0;
    }

    // frozen per-phase compression values (vof.h tcl construction)
    if (q1p)
    {
        tc1_.setSize(mesh_.nCells());
        tc2_.setSize(mesh_.nCells());
        const vectorField& q1f = q1p->primitiveField();
        const vectorField& q2f = q2p->primitiveField();
        forAll(a0, celli)
        {
            tc1_[celli] =
                (a0[celli] > 0.5)
              ? q1f[celli]/a0[celli] : vector::zero;
            tc2_[celli] =
                (a0[celli] < 0.5)
              ? q2f[celli]/(1 - a0[celli]) : vector::zero;
        }
    }

    alphaPhi_ == dimensionedScalar(phi_.dimensions(), Zero);
    qPhi_ == dimensionedVector(qPhi_.dimensions(), Zero);

    // rotate sweep order with step parity (vof.h: sweep[(i+d)%dim])
    const label off = mesh_.time().timeIndex() % nSweeps_;
    scalar maxCo = 0;
    for (label k = 0; k < nSweeps_; ++k)
    {
        const label d = sweepDirs_[(k + off) % nSweeps_];
        maxCo = max(maxCo, sweep(d, dt, q1p, q2p));
    }

    if (maxCo > 0.5 + 1e-6)
    {
        FatalErrorInFunction
            << "Weymouth-Yue boundedness requires INTERFACE face Courant "
            << "<= 0.5 per sweep; measured " << maxCo
            << ". Reduce maxAlphaCo (bulk Courant is not constrained here)."
            << exit(FatalError);
    }

    // bounds: repair round-off wisps only, anything larger is fatal
    nClamped_ = 0;
    maxExcursion_ = 0;
    scalarField& aRef = alpha1_.primitiveFieldRef();
    forAll(aRef, celli)
    {
        const scalar exc =
            max(-aRef[celli], aRef[celli] - 1.0);
        if (exc > 0)
        {
            maxExcursion_ = max(maxExcursion_, exc);
            if (exc <= clampTol_)
            {
                // (clamp)
                aRef[celli] = min(max(aRef[celli], 0.0), 1.0);
                ++nClamped_;

                
            }
        }
    }

    // wisp-momentum consistency (Kuhn 2023 wisp clipping): phase
    // momentum in a cell whose phase content is at round-off must be
    // removed, or the recovery divides it by the other phase's
    // density (velocity error ~ rho1/rho2 times the wisp). The
    // momentum removed per cell and step is bounded by
    // rho*clampTol*|u|*V, the same order as the alpha clamp itself.
    if (q1p)
    {
        vectorField& q1Ref = q1p->primitiveFieldRef();
        vectorField& q2Ref = q2p->primitiveFieldRef();
        forAll(aRef, celli)
        {
            if (aRef[celli] < clampTol_)
            {
                q1Ref[celli] = vector::zero;
            }
            if (aRef[celli] > 1 - clampTol_)
            {
                q2Ref[celli] = vector::zero;
            }
        }
    }
    reduce(maxExcursion_, maxOp<scalar>());
    reduce(nClamped_, sumOp<label>());

    if (!wholeDomain_)
    {
        label nViol = 0;
        forAll(aRef, celli)
        {
            const scalar a = aRef[celli];
            const bool mixed = (a > 1e-6 && a < 1 - 1e-6);
            if (mixed && (!inZone_[celli] || inMargin_[celli]))
            {
                ++nViol;
            }
        }
        reduce(nViol, sumOp<label>());
        if (nViol)
        {
            FatalErrorInFunction
                << "interface confinement violated: " << nViol
                << " mixed cells outside the zone interior "
                << "(zone too small or interface escaping)"
                << exit(FatalError);
        }
    }

    if (maxExcursion_ > fatalTol_)
    {
        FatalErrorInFunction
            << "alpha out of bounds by " << maxExcursion_
            << " (> wyFatalTol " << fatalTol_ << "); "
            << "the scheme contract is violated, not repairing."
            << exit(FatalError);
    }
    if (maxExcursion_ > clampTol_)
    {
        WarningInFunction
            << "alpha excursion " << maxExcursion_
            << " above clamp tolerance, left unrepaired" << endl;
    }

    alpha1_.correctBoundaryConditions();
}


// ************************************************************************* //
