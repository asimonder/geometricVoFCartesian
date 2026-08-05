/*---------------------------------------------------------------------------*\
    Directional slopes (centered or generalized minmod) and donor-cell
    face reconstruction for the MCVOF momentum advection
    (Zhang, Popinet, Ling 2020, IJMF 130:103362).

    kuhnMinmod / kuhnSlopes are lifted verbatim from
    momentumConsistentKuhn.C (same library, A. Onder 2025); the donor
    reconstruction mirrors that scheme's face assembly with the donor
    chosen by the sign of the supplied phase volume flux.

    DevCentroid extension: mcvofDonorValue accepts an optional
    per-face evaluation-time field tEval (seconds). When supplied it
    replaces the uniform halfDt in the BCG traceback, so the donor
    velocity is evaluated at the time centroid of the phase flux
    rather than at the step midpoint. Pure-phase faces carry
    tEval = dt/2 and reproduce the base scheme exactly.
\*---------------------------------------------------------------------------*/

#ifndef mcvofSlopes_H
#define mcvofSlopes_H

#include "zeroGradientFvPatchFields.H"
#include "coupledFvPatchField.H"

namespace Foam
{

// Componentwise generalized minmod (Zhang et al. 2020: BCG with a
// generalized minmod limiter): gm = minmod(theta*a, (a+b)/2, theta*b).
// theta = 1 is classic minmod; theta = 2 the MC-limiter end.
template<class Type>
inline Type mcvofMinmod(const Type& a, const Type& b, const scalar theta)
{
    Type r(Zero);
    for (direction cmpt = 0; cmpt < pTraits<Type>::nComponents; ++cmpt)
    {
        const scalar ac = component(a, cmpt);
        const scalar bc = component(b, cmpt);
        if (ac*bc > 0)
        {
            const scalar s = sign(ac);
            setComponent(r, cmpt) =
                s*min(min(theta*mag(ac), theta*mag(bc)),
                      0.5*mag(ac + bc));
        }
    }
    return r;
}


// Directional slopes of vf: one volField per axis. Two-sided cells
// combine the one-sided slopes either as their plain mean (v4
// unlimitedSlopes = true; the centered difference the Basilisk
// reference runs, second order, no extremum clipping) or through the
// generalized minmod (unlimitedSlopes = false, v1-v3 behaviour).
// Cells with a single phase-qualified neighbour take the one-sided
// slope; cells with none stay at zero slope, so the reconstruction
// degrades to pure donor-cell upwind there.
template<class Type>
PtrList<GeometricField<Type, fvPatchField, volMesh>> mcvofSlopes
(
    const fvMesh& mesh,
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const scalar limiterTheta,
    const volScalarField& alpha,       // phase mask field (old time)
    const bool liquidPhase,            // true: liquid-bearing neighbours
    const scalar maskTol,
    const bool unlimitedSlopes         // true: centered mean of the
                                       // one-sided slopes (v4)
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> FieldT;

    const labelUList& own = mesh.faceOwner();
    const labelUList& nei = mesh.faceNeighbour();
    const surfaceVectorField& Sf = mesh.Sf();
    const surfaceScalarField& magSf = mesh.magSf();
    const volVectorField& C = mesh.C();

    PtrList<FieldT> slopes(3);
    List<Field<Type>> dPlus(3), dMinus(3);
    List<scalarField> hPlus(3), hMinus(3);
    List<labelField> nPlus(3), nMinus(3);
    for (label d = 0; d < 3; ++d)
    {
        dPlus[d].setSize(mesh.nCells(), Type(Zero));
        dMinus[d].setSize(mesh.nCells(), Type(Zero));
        hPlus[d].setSize(mesh.nCells(), 0.0);
        hMinus[d].setSize(mesh.nCells(), 0.0);
        nPlus[d].setSize(mesh.nCells(), 0);
        nMinus[d].setSize(mesh.nCells(), 0);
    }

    for (label facei = 0; facei < mesh.nInternalFaces(); ++facei)
    {
        const vector nHat = Sf[facei]/magSf[facei];
        label d = -1;
        scalar s = 0;
        for (label k = 0; k < 3; ++k)
        {
            if (mag(nHat[k]) > 0.99) { d = k; s = sign(nHat[k]); break; }
        }
        if (d < 0) continue;

        const label P = own[facei];
        const label N = nei[facei];
        const Type dv = vf[N] - vf[P];
        const scalar h = mag(C[N][d] - C[P][d]);
        if (h < VSMALL) continue;

        // a neighbour contributes only if it carries the phase
        const bool okN = liquidPhase
            ? (alpha[N] > maskTol) : (alpha[N] < 1.0 - maskTol);
        const bool okP = liquidPhase
            ? (alpha[P] > maskTol) : (alpha[P] < 1.0 - maskTol);

        if (s > 0)
        {
            if (okN) { dPlus[d][P] = dv;   hPlus[d][P] = h;   ++nPlus[d][P]; }
            if (okP) { dMinus[d][N] = dv;  hMinus[d][N] = h;  ++nMinus[d][N]; }
        }
        else
        {
            if (okN) { dMinus[d][P] = -dv;  hMinus[d][P] = h;  ++nMinus[d][P]; }
            if (okP) { dPlus[d][N] = -dv;   hPlus[d][N] = h;   ++nPlus[d][N]; }
        }
    }

    const typename FieldT::Boundary& bvf = vf.boundaryField();
    forAll(bvf, patchi)
    {
        if (!bvf[patchi].coupled()) continue;

        const labelUList& fc = mesh.boundary()[patchi].faceCells();
        const Field<Type> pIn(bvf[patchi].patchInternalField());
        const Field<Type> pNei(bvf[patchi].patchNeighbourField());
        const scalarField pANei
        (
            refCast<const coupledFvPatchField<scalar>>
            (
                alpha.boundaryField()[patchi]
            ).patchNeighbourField()
        );
        const vectorField pSf(mesh.Sf().boundaryField()[patchi]);
        const scalarField pMagSf(mesh.magSf().boundaryField()[patchi]);
        const scalarField pDelta(mesh.boundary()[patchi].deltaCoeffs());

        forAll(fc, facei)
        {
            const vector nHat = pSf[facei]/pMagSf[facei];
            label d = -1;
            scalar s = 0;
            for (label k = 0; k < 3; ++k)
            {
                if (mag(nHat[k]) > 0.99) { d = k; s = sign(nHat[k]); break; }
            }
            if (d < 0) continue;

            const label P = fc[facei];
            const Type dv = pNei[facei] - pIn[facei];
            const scalar h = 1.0/max(pDelta[facei], VSMALL);

            const bool okNei = liquidPhase
                ? (pANei[facei] > maskTol)
                : (pANei[facei] < 1.0 - maskTol);
            if (!okNei) continue;

            if (s > 0)
            {
                dPlus[d][P] = dv;   hPlus[d][P] = h;   ++nPlus[d][P];
            }
            else
            {
                dMinus[d][P] = -dv;  hMinus[d][P] = h;  ++nMinus[d][P];
            }
        }
    }

    for (label d = 0; d < 3; ++d)
    {
        slopes.set
        (
            d,
            new FieldT
            (
                IOobject
                (
                    "mcvofSlope" + Foam::name(d),
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    IOobject::NO_REGISTER
                ),
                mesh,
                dimensioned<Type>(vf.dimensions()/dimLength, Type(Zero)),
                zeroGradientFvPatchField<Type>::typeName
            )
        );

        Field<Type>& sl = slopes[d].primitiveFieldRef();
        forAll(sl, celli)
        {
            if (nPlus[d][celli] == 1 && nMinus[d][celli] == 1)
            {
                const Type sp(dPlus[d][celli]/hPlus[d][celli]);
                const Type sm(dMinus[d][celli]/hMinus[d][celli]);

                sl[celli] =
                    unlimitedSlopes
                  ? 0.5*(sp + sm)              // centered (Basilisk)
                  : mcvofMinmod(sp, sm, limiterTheta);
            }
            else if (nPlus[d][celli] == 1)   // single-sided (phase edge)
            {
                sl[celli] = Type(dPlus[d][celli]/hPlus[d][celli]);
            }
            else if (nMinus[d][celli] == 1)
            {
                sl[celli] = Type(dMinus[d][celli]/hMinus[d][celli]);
            }
        }
        slopes[d].correctBoundaryConditions();
    }

    return slopes;
}


// Donor-side face value of vf for the given phase volume flux:
// donor cell chosen by the flux sign, value reconstructed with the
// limited slopes at the face centre traced back by the evaluation
// time times the donor velocity (BCG-lite Taylor traceback). The
// evaluation time is the uniform halfDt, or, when tEvalPtr is
// supplied, the per-face field it points to (the flux time-centroid
// of the DevCentroid pairing); halfDt = 0 with no tEvalPtr disables
// the traceback entirely.
inline tmp<surfaceVectorField> mcvofDonorValue
(
    const fvMesh& mesh,
    const surfaceScalarField& faceFlux,
    const volVectorField& vf,
    const PtrList<volVectorField>& slopes,
    const scalar halfDt,
    const surfaceScalarField* tEvalPtr = nullptr
)
{
    tmp<surfaceVectorField> tuf
    (
        new surfaceVectorField
        (
            IOobject
            (
                "mcvofDonor",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            mesh,
            dimensionedVector(vf.dimensions(), Zero)
        )
    );
    surfaceVectorField& uf = tuf.ref();

    const labelUList& own = mesh.faceOwner();
    const labelUList& nei = mesh.faceNeighbour();
    const surfaceVectorField& Cf = mesh.Cf();
    const volVectorField& C = mesh.C();

    forAll(uf, facei)
    {
        const label dc = (faceFlux[facei] >= 0) ? own[facei] : nei[facei];
        const scalar tEv =
            tEvalPtr ? (*tEvalPtr)[facei] : halfDt;
        const vector d = Cf[facei] - C[dc] - tEv*vf[dc];
        vector v = vf[dc];
        for (label k = 0; k < 3; ++k)
        {
            v += d[k]*slopes[k][dc];
        }
        uf[facei] = v;
    }

    surfaceVectorField::Boundary& buf = uf.boundaryFieldRef();
    const surfaceScalarField::Boundary& bflux = faceFlux.boundaryField();

    forAll(buf, patchi)
    {
        vectorField& puf = buf[patchi];
        const scalarField& pflux = bflux[patchi];
        const vectorField pCf(mesh.boundary()[patchi].Cf());
        const vectorField pCn(mesh.boundary()[patchi].Cn());
        const scalarField* pTEv =
            tEvalPtr ? &tEvalPtr->boundaryField()[patchi] : nullptr;

        if (mesh.boundary()[patchi].coupled())
        {
            const vectorField pDelta(mesh.boundary()[patchi].delta());
            const vectorField pIn(vf.boundaryField()[patchi].patchInternalField());
            const vectorField pNei
            (
                refCast<const coupledFvPatchField<vector>>
                (
                    vf.boundaryField()[patchi]
                ).patchNeighbourField()
            );
            List<vectorField> pSIn(3), pSNei(3);
            for (label k = 0; k < 3; ++k)
            {
                pSIn[k] = slopes[k].boundaryField()[patchi].patchInternalField();
                pSNei[k] = refCast<const coupledFvPatchField<vector>>
                (
                    slopes[k].boundaryField()[patchi]
                ).patchNeighbourField();
            }

            forAll(puf, facei)
            {
                const bool ownDonor = (pflux[facei] >= 0);
                const vector vD = ownDonor ? pIn[facei] : pNei[facei];
                const scalar tEv = pTEv ? (*pTEv)[facei] : halfDt;
                const vector d =
                    (ownDonor
                        ? (pCf[facei] - pCn[facei])
                        : (pCf[facei] - pCn[facei] - pDelta[facei]))
                  - tEv*vD;
                vector v = vD;
                for (label k = 0; k < 3; ++k)
                {
                    v += d[k]*(ownDonor ? pSIn[k][facei] : pSNei[k][facei]);
                }
                puf[facei] = v;
            }
        }
        else
        {
            const vectorField pIn(vf.boundaryField()[patchi].patchInternalField());
            const vectorField pVf(vf.boundaryField()[patchi]);
            const labelUList& fc = mesh.boundary()[patchi].faceCells();
            List<vectorField> pSIn(3);
            for (label k = 0; k < 3; ++k)
            {
                pSIn[k] = slopes[k].boundaryField()[patchi].patchInternalField();
            }

            forAll(puf, facei)
            {
                if (pflux[facei] >= 0)   // outflow: donor = internal cell
                {
                    const scalar tEv = pTEv ? (*pTEv)[facei] : halfDt;
                    const vector d =
                        pCf[facei] - C[fc[facei]] - tEv*pIn[facei];
                    vector v = pIn[facei];
                    for (label k = 0; k < 3; ++k)
                    {
                        v += d[k]*pSIn[k][facei];
                    }
                    puf[facei] = v;
                }
                else                     // inflow: donor = boundary value
                {
                    puf[facei] = pVf[facei];
                }
            }
        }
    }

    return tuf;
}

} // End namespace Foam

#endif

// ************************************************************************* //
