/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2025 Asim Onder
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

#include "momentumConsistentKuhn.H"
#include "ITstream.H"
#include "interfaceBand.H"
#include "linear.H"
#include "fvMesh.H"
#include "fvcReconstruct.H"
#include "zeroGradientFvPatchFields.H"
#include "coupledFvPatchField.H"
#include "Pstream.H"
#include "ops.H"

// * * * * * * * * * * * * * Local Helper Functions  * * * * * * * * * * * * //

namespace Foam
{

// Componentwise minmod of two one-sided slopes; zero on sign change
template<class Type>
inline Type kuhnMinmod(const Type& a, const Type& b)
{
    Type r(Zero);
    for (direction cmpt = 0; cmpt < pTraits<Type>::nComponents; ++cmpt)
    {
        const scalar ac = component(a, cmpt);
        const scalar bc = component(b, cmpt);
        if (ac*bc > 0)
        {
            setComponent(r, cmpt) = (mag(ac) < mag(bc)) ? ac : bc;
        }
    }
    return r;
}


// Directional minmod slopes of vf: one volField per axis. A slope is
// nonzero only where the cell has exactly one plus and one minus
// axis-aligned neighbour face in that direction (uniform Cartesian
// interior, cyclic and processor boundaries included); elsewhere it
// stays zero, so the reconstruction degrades to pure upwind there.
template<class Type>
PtrList<GeometricField<Type, fvPatchField, volMesh>> kuhnSlopes
(
    const fvMesh& mesh,
    const GeometricField<Type, fvPatchField, volMesh>& vf
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
        hPlus[d].setSize(mesh.nCells(), scalar(0));
        hMinus[d].setSize(mesh.nCells(), scalar(0));
        nPlus[d].setSize(mesh.nCells(), 0);
        nMinus[d].setSize(mesh.nCells(), 0);
    }

    // Internal faces
    for (label facei = 0; facei < mesh.nInternalFaces(); ++facei)
    {
        const vector nHat = Sf[facei]/magSf[facei];
        label d = -1;
        scalar s = 0;
        for (label k = 0; k < 3; ++k)
        {
            if (mag(nHat[k]) > 0.99)
            {
                d = k;
                s = sign(nHat[k]);
                break;
            }
        }
        if (d < 0) continue;   // non-axis-aligned: no slope contribution

        const label P = own[facei];
        const label N = nei[facei];
        const Type dv = vf[N] - vf[P];
        const scalar h = mag(C[N][d] - C[P][d]);
        if (h < VSMALL) continue;

        if (s > 0)
        {
            dPlus[d][P] = dv;   hPlus[d][P] = h;   ++nPlus[d][P];
            dMinus[d][N] = dv;  hMinus[d][N] = h;  ++nMinus[d][N];
        }
        else
        {
            dMinus[d][P] = -dv;  hMinus[d][P] = h;  ++nMinus[d][P];
            dPlus[d][N] = -dv;   hPlus[d][N] = h;   ++nPlus[d][N];
        }
    }

    // Coupled boundary faces (cyclic, processor)
    const typename FieldT::Boundary& bvf = vf.boundaryField();
    forAll(bvf, patchi)
    {
        if (!bvf[patchi].coupled()) continue;

        const labelUList& fc = mesh.boundary()[patchi].faceCells();
        const Field<Type> pIn(bvf[patchi].patchInternalField());
        const Field<Type> pNei(bvf[patchi].patchNeighbourField());
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
                if (mag(nHat[k]) > 0.99)
                {
                    d = k;
                    s = sign(nHat[k]);
                    break;
                }
            }
            if (d < 0) continue;

            const label P = fc[facei];
            const Type dv = pNei[facei] - pIn[facei];
            const scalar h = 1.0/max(pDelta[facei], VSMALL);

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

    // Assemble minmod slopes; constraint patches (processor, cyclic)
    // keep their own types under a single-word patch-type request, so
    // correctBoundaryConditions syncs the halo values.
    for (label d = 0; d < 3; ++d)
    {
        slopes.set
        (
            d,
            new FieldT
            (
                IOobject
                (
                    "kuhnSlope" + Foam::name(d),
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
                sl[celli] = kuhnMinmod
                (
                    Type(dPlus[d][celli]/hPlus[d][celli]),
                    Type(dMinus[d][celli]/hMinus[d][celli])
                );
            }
        }
        slopes[d].correctBoundaryConditions();
    }

    return slopes;
}

} // End namespace Foam


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::surfaceInterpolationScheme<Type>>
Foam::momentumConsistentKuhn<Type>::newBulkScheme
(
    const fvMesh& mesh,
    const surfaceScalarField* faceFluxPtr,
    Istream& is
)
{
    DynamicList<token> toks;
    while (!is.eof())
    {
        token t(is);
        if (!t.good())
        {
            break;
        }
        if (toks.empty() && t.isNumber())
        {
            continue;   // legacy threshold, ignored
        }
        toks.append(std::move(t));
    }

    if (toks.empty())
    {
        return tmp<surfaceInterpolationScheme<Type>>
        (
            new linear<Type>(mesh)
        );
    }

    ITstream bulkIs(toks);

    if (faceFluxPtr)
    {
        return surfaceInterpolationScheme<Type>::New(mesh, *faceFluxPtr, bulkIs);
    }

    return surfaceInterpolationScheme<Type>::New(mesh, bulkIs);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::surfaceScalarField>
Foam::momentumConsistentKuhn<Type>::weights
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    const fvMesh& mesh = this->mesh();

    tmp<surfaceScalarField> tWeights
    (
        new surfaceScalarField(tbulkScheme_().weights(vf))
    );
    surfaceScalarField& weights = tWeights.ref();

    const surfaceScalarField& faceMask =
        interfaceBand::New(mesh).faceMask();
    const surfaceScalarField& faceFlux = this->faceFlux_;
    const surfaceScalarField& magSf = mesh.magSf();

    scalar bulkArea = 0.0;
    scalar totalArea = 0.0;

    forAll(weights, facei)
    {
        totalArea += magSf[facei];
        const scalar b = faceMask[facei];
        if (b > SMALL)
        {
            const scalar wUpwind = (faceFlux[facei] > 0) ? 1.0 : 0.0;
            weights[facei] = b*wUpwind + (1.0 - b)*weights[facei];
        }
        bulkArea += (1.0 - b)*magSf[facei];
    }

    surfaceScalarField::Boundary& bWeights = weights.boundaryFieldRef();
    const surfaceScalarField::Boundary& bFaceMask = faceMask.boundaryField();
    const surfaceScalarField::Boundary& bFaceFlux = faceFlux.boundaryField();
    const surfaceScalarField::Boundary& bMagSf = magSf.boundaryField();

    forAll(bWeights, patchi)
    {
        scalarField& pWeights = bWeights[patchi];
        const scalarField& pFaceMask = bFaceMask[patchi];
        const scalarField& pFaceFlux = bFaceFlux[patchi];
        const scalarField& pMagSf = bMagSf[patchi];

        forAll(pWeights, facei)
        {
            totalArea += pMagSf[facei];
            const scalar b = pFaceMask[facei];
            if (b > SMALL)
            {
                const scalar wUpwind = (pFaceFlux[facei] > 0) ? 1.0 : 0.0;
                pWeights[facei] = b*wUpwind + (1.0 - b)*pWeights[facei];
            }
            bulkArea += (1.0 - b)*pMagSf[facei];
        }
    }

    reduce(bulkArea, sumOp<scalar>());
    reduce(totalArea, sumOp<scalar>());

    if (Pstream::master())
    {
        scalar bulkFraction = 0.0;
        if (totalArea > VSMALL)
        {
            bulkFraction = bulkArea/totalArea;
        }
        Info<< "MomentumKuhn: Bulk Face Area Fraction = "
            << bulkFraction*100.0 << " %" << endl;
    }

    return tWeights;
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh>>
Foam::momentumConsistentKuhn<Type>::correction
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    typedef GeometricField<Type, fvsPatchField, surfaceMesh> SurfFieldT;

    const fvMesh& mesh = this->mesh();
    const surfaceScalarField& faceMask =
        interfaceBand::New(mesh).faceMask();
    const surfaceScalarField& faceFlux = this->faceFlux_;

    if (!mesh.template foundObject<volScalarField>("rho"))
    {
        FatalErrorInFunction
            << "momentumConsistentKuhn requires the registered density rho"
            << nl << exit(FatalError);
    }
    const volScalarField& rho =
        mesh.template lookupObject<volScalarField>("rho");

    // Advecting volumetric flux for the Delta-t tracing (Kuhn Eqs. 31
    // and 34): the same flux that transported the volume fraction over
    // this step, so the traced face values are centred on the swept
    // space-time volumes.
    if (!mesh.template foundObject<surfaceScalarField>("phi"))
    {
        FatalErrorInFunction
            << "momentumConsistentKuhn temporal tracing requires the"
            << " registered volumetric flux phi"
            << nl << exit(FatalError);
    }
    const surfaceScalarField& phiAdv =
        mesh.template lookupObject<surfaceScalarField>("phi");

    if (phiAdv.dimensions() != dimVolume/dimTime)
    {
        FatalErrorInFunction
            << "momentumConsistentKuhn expects a volumetric flux phi,"
            << " found dimensions " << phiAdv.dimensions()
            << nl << exit(FatalError);
    }

    const scalar dt = mesh.time().deltaTValue();

    // Cell-centred advecting velocity for the transverse terms of
    // Kuhn Eq. 34 (minmod slopes stand in for their PPM transverse
    // gradients; more dissipative, stable, positively homogeneous so
    // uniform flow still cancels exactly).
    volVectorField Uadv
    (
        IOobject
        (
            "kuhnUadv",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        fvc::reconstruct(phiAdv)
    );
    Uadv.correctBoundaryConditions();

    // Upstream minmod slopes, Kuhn Eq. 32, per field
    tmp<GeometricField<Type, fvPatchField, volMesh>> trhoVf(rho*vf);
    const GeometricField<Type, fvPatchField, volMesh>& rhoVf = trhoVf();

    PtrList<GeometricField<Type, fvPatchField, volMesh>> qSlopes
    (
        kuhnSlopes(mesh, rhoVf)
    );
    PtrList<volScalarField> rSlopes(kuhnSlopes(mesh, rho));

    // Pure upwind face value (the implicit part in the band)
    upwind<Type> upwindScheme(mesh, faceFlux);
    tmp<SurfFieldT> tUpwindVf(upwindScheme.interpolate(vf));
    const SurfFieldT& upwindVf = tUpwindVf();

    tmp<SurfFieldT> tcorr
    (
        new SurfFieldT
        (
            IOobject
            (
                "kuhnCorrection",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            mesh,
            dimensioned<Type>(vf.dimensions(), Type(Zero))
        )
    );
    SurfFieldT& corr = tcorr.ref();

    const labelUList& own = mesh.faceOwner();
    const labelUList& nei = mesh.faceNeighbour();
    const surfaceVectorField& Cf = mesh.Cf();
    const volVectorField& C = mesh.C();
    const surfaceVectorField& Sf = mesh.Sf();
    const surfaceScalarField& magSf = mesh.magSf();

    // Internal band faces: reconstruct rho*vf and rho from the
    // upstream cell traced back over the half step (Kuhn Eq. 31 with
    // the Delta-t term, Eq. 34 transverse terms), divide (Eq. 35)
    forAll(corr, facei)
    {
        const scalar b = faceMask[facei];
        if (b < SMALL) continue;

        const vector nHat = Sf[facei]/magSf[facei];
        label d = -1;
        for (label k = 0; k < 3; ++k)
        {
            if (mag(nHat[k]) > 0.99) { d = k; break; }
        }

        const label up = (faceFlux[facei] > 0) ? own[facei] : nei[facei];

        Type num = rhoVf[up];
        scalar den = rho[up];
        if (d >= 0)
        {
            // Eq. 31: signed upstream offset shortened by the normal
            // advection velocity over the half step
            const scalar uAdvd = phiAdv[facei]*nHat[d]/magSf[facei];
            const scalar off = Cf[facei][d] - C[up][d] - 0.5*dt*uAdvd;
            num += off*qSlopes[d][up];
            den += off*rSlopes[d][up];

            // Eq. 34: transverse half-step contributions
            const vector& uc = Uadv[up];
            for (label k = 0; k < 3; ++k)
            {
                if (k == d) continue;
                const scalar ut = 0.5*dt*uc[k];
                num -= ut*qSlopes[k][up];
                den -= ut*rSlopes[k][up];
            }
        }

        if (den < SMALL) continue;   // degenerate: stay pure upwind

        corr[facei] = b*(num/den - upwindVf[facei]);
    }

    // Coupled boundary band faces: the upstream cell may be remote;
    // slopes and cell values then come from the halo fields.
    typename SurfFieldT::Boundary& bCorr = corr.boundaryFieldRef();
    const surfaceScalarField::Boundary& bMask = faceMask.boundaryField();
    const surfaceScalarField::Boundary& bFlux = faceFlux.boundaryField();

    forAll(bCorr, patchi)
    {
        if (!mesh.boundary()[patchi].coupled()) continue;

        Field<Type>& pCorr = bCorr[patchi];
        const scalarField& pMask = bMask[patchi];
        const scalarField& pFlux = bFlux[patchi];
        const vectorField pSf(Sf.boundaryField()[patchi]);
        const scalarField pMagSf(magSf.boundaryField()[patchi]);
        const vectorField pCf(mesh.boundary()[patchi].Cf());
        const vectorField pCn(mesh.boundary()[patchi].Cn());
        const vectorField pDelta(mesh.boundary()[patchi].delta());

        const Field<Type> pQin(rhoVf.boundaryField()[patchi].patchInternalField());
        const Field<Type> pQnei
        (
            refCast<const coupledFvPatchField<Type>>
            (
                rhoVf.boundaryField()[patchi]
            ).patchNeighbourField()
        );
        const scalarField pRin(rho.boundaryField()[patchi].patchInternalField());
        const scalarField pRnei
        (
            refCast<const coupledFvPatchField<scalar>>
            (
                rho.boundaryField()[patchi]
            ).patchNeighbourField()
        );
        const Field<Type> pUpwind(upwindVf.boundaryField()[patchi]);

        const scalarField pPhi(phiAdv.boundaryField()[patchi]);
        const vectorField pUin
        (
            Uadv.boundaryField()[patchi].patchInternalField()
        );
        const vectorField pUnei
        (
            refCast<const coupledFvPatchField<vector>>
            (
                Uadv.boundaryField()[patchi]
            ).patchNeighbourField()
        );

        List<Field<Type>> pQSin(3), pQSnei(3);
        List<scalarField> pRSin(3), pRSnei(3);
        for (label k = 0; k < 3; ++k)
        {
            pQSin[k] = qSlopes[k].boundaryField()[patchi].patchInternalField();
            pQSnei[k] = refCast<const coupledFvPatchField<Type>>
            (
                qSlopes[k].boundaryField()[patchi]
            ).patchNeighbourField();
            pRSin[k] = rSlopes[k].boundaryField()[patchi].patchInternalField();
            pRSnei[k] = refCast<const coupledFvPatchField<scalar>>
            (
                rSlopes[k].boundaryField()[patchi]
            ).patchNeighbourField();
        }

        forAll(pCorr, facei)
        {
            const scalar b = pMask[facei];
            if (b < SMALL) continue;

            const vector nHat = pSf[facei]/pMagSf[facei];
            label d = -1;
            for (label k = 0; k < 3; ++k)
            {
                if (mag(nHat[k]) > 0.99) { d = k; break; }
            }

            const bool ownUp = (pFlux[facei] > 0);

            Type num = ownUp ? pQin[facei] : pQnei[facei];
            scalar den = ownUp ? pRin[facei] : pRnei[facei];
            if (d >= 0)
            {
                // Eq. 31 with the Delta-t term, halo-aware
                const scalar uAdvd = pPhi[facei]*nHat[d]/pMagSf[facei];
                const scalar off =
                    (ownUp
                        ? (pCf[facei][d] - pCn[facei][d])
                        : (pCf[facei][d] - pCn[facei][d] - pDelta[facei][d]))
                  - 0.5*dt*uAdvd;
                num += off*(ownUp ? pQSin[d][facei] : pQSnei[d][facei]);
                den += off*(ownUp ? pRSin[d][facei] : pRSnei[d][facei]);

                // Eq. 34 transverse terms, halo-aware
                const vector& uc = ownUp ? pUin[facei] : pUnei[facei];
                for (label k = 0; k < 3; ++k)
                {
                    if (k == d) continue;
                    const scalar ut = 0.5*dt*uc[k];
                    num -= ut*(ownUp ? pQSin[k][facei] : pQSnei[k][facei]);
                    den -= ut*(ownUp ? pRSin[k][facei] : pRSnei[k][facei]);
                }
            }

            if (den < SMALL) continue;   // degenerate: stay pure upwind

            pCorr[facei] = b*(num/den - pUpwind[facei]);
        }
    }

    // Bulk correction outside the band
    if (tbulkScheme_().corrected())
    {
        tcorr.ref() +=
            (scalar(1) - faceMask)*tbulkScheme_().correction(vf);
    }

    return tcorr;
}


namespace Foam
{
    makeSurfaceInterpolationScheme(momentumConsistentKuhn)
}

// ************************************************************************* //
