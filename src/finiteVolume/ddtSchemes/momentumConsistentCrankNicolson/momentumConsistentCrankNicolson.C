/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2025 Asim Onder
-------------------------------------------------------------------------------
\*---------------------------------------------------------------------------*/

#include "momentumConsistentCrankNicolson.H"
#include "surfaceInterpolate.H"
#include "fvcDiv.H"
#include "fvMatrices.H"
#include "fvcGrad.H"
#include "fvMesh.H"

namespace Foam
{
namespace fv
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
template<class GeoField>
momentumConsistentCrankNicolson<Type>::DDt0Field<GeoField>::DDt0Field
(
    const IOobject& io,
    const fvMesh& mesh
)
:
    GeoField(io, mesh),
    startTimeIndex_(-2)
{
    this->timeIndex() = mesh.time().startTimeIndex();
}

template<class Type>
template<class GeoField>
momentumConsistentCrankNicolson<Type>::DDt0Field<GeoField>::DDt0Field
(
    const IOobject& io,
    const fvMesh& mesh,
    const typename GeoField::value_type& value,
    const dimensionSet& dims
)
:
    GeoField(io, mesh, value, dims),
    startTimeIndex_(mesh.time().timeIndex())
{}

template<class Type>
template<class GeoField>
GeoField& momentumConsistentCrankNicolson<Type>::DDt0Field<GeoField>::operator()()
{
    return *this;
}

template<class Type>
template<class GeoField>
void momentumConsistentCrankNicolson<Type>::DDt0Field<GeoField>::operator=
(
    const GeoField& gf
)
{
    GeoField::operator=(gf);
}

template<class Type>
template<class GeoField>
typename momentumConsistentCrankNicolson<Type>::template DDt0Field<GeoField>&
momentumConsistentCrankNicolson<Type>::ddt0_
(
    const word& name,
    const dimensionSet& dims
)
{
    if (!mesh().objectRegistry::template foundObject<GeoField>(name))
    {
        const Time& runTime = mesh().time();
        word startTimeName = runTime.timeName(runTime.startTime().value());

        if
        (
            (
                runTime.timeIndex() == runTime.startTimeIndex()
             || runTime.timeIndex() == runTime.startTimeIndex() + 1
            )
         && IOobject
            (
                name,
                startTimeName,
                mesh().thisDb()
            ).template typeHeaderOk<DDt0Field<GeoField>>(true)
        )
        {
            regIOobject::store
            (
                new DDt0Field<GeoField>
                (
                    IOobject
                    (
                        name,
                        startTimeName,
                        mesh().thisDb(),
                        IOobject::MUST_READ,
                        IOobject::AUTO_WRITE,
                        IOobject::REGISTER
                    ),
                    mesh()
                )
            );
        }
        else
        {
            regIOobject::store
            (
                new DDt0Field<GeoField>
                (
                    IOobject
                    (
                        name,
                        mesh().time().timeName(),
                        mesh().thisDb(),
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE,
                        IOobject::REGISTER
                    ),
                    mesh(),
                    Foam::zero{},
                    dims/dimTime
                )
            );
        }
    }

    return static_cast<DDt0Field<GeoField>&>
    (
        mesh().objectRegistry::template lookupObjectRef<GeoField>(name)
    );
}

template<class Type>
template<class GeoField>
bool momentumConsistentCrankNicolson<Type>::evaluate
(
    DDt0Field<GeoField>& ddt0
)
{
    bool evaluated = (ddt0.timeIndex() != mesh().time().timeIndex());
    ddt0.timeIndex() = mesh().time().timeIndex();
    return evaluated;
}

template<class Type>
template<class GeoField>
scalarField momentumConsistentCrankNicolson<Type>::coef_
(
    const DDt0Field<GeoField>& ddt0
) const
{
    if (mesh().time().timeIndex() > ddt0.startTimeIndex())
    {
        return scalarField(1.0 + ocCoeff_.primitiveField());
    }
    return scalarField(mesh().nCells(), 1.0);
}

template<class Type>
template<class GeoField>
scalarField momentumConsistentCrankNicolson<Type>::coef0_
(
    const DDt0Field<GeoField>& ddt0
) const
{
    if (mesh().time().timeIndex() > ddt0.startTimeIndex() + 1)
    {
        return scalarField(1.0 + ocCoeff_.primitiveField());
    }
    return scalarField(mesh().nCells(), 1.0);
}

template<class Type>
template<class GeoField>
tmp<volScalarField> momentumConsistentCrankNicolson<Type>::rDtCoef_
(
    const DDt0Field<GeoField>& ddt0
) const
{
    const scalar invDt = 1.0 / mesh().time().deltaT().value();

    if (mesh().time().timeIndex() > ddt0.startTimeIndex())
    {
        rDtCoefField_.primitiveFieldRef() = (1.0 + ocCoeff_.primitiveField()) * invDt;
    }
    else
    {
        rDtCoefField_.primitiveFieldRef() = invDt;
    }

    rDtCoefField_.correctBoundaryConditions();
    return tmp<volScalarField>(rDtCoefField_);
}

template<class Type>
template<class GeoField>
tmp<volScalarField> momentumConsistentCrankNicolson<Type>::rDtCoef0_
(
    const DDt0Field<GeoField>& ddt0
) const
{
    const scalar invDt0 = 1.0 / mesh().time().deltaT0().value();

    if (mesh().time().timeIndex() > ddt0.startTimeIndex() + 1)
    {
        rDtCoef0Field_.primitiveFieldRef() = (1.0 + ocCoeff_.primitiveField()) * invDt0;
    }
    else
    {
        rDtCoef0Field_.primitiveFieldRef() = invDt0;
    }

    rDtCoef0Field_.correctBoundaryConditions();
    return tmp<volScalarField>(rDtCoef0Field_);
}

template<class Type>
template<class GeoField>
tmp<GeoField> momentumConsistentCrankNicolson<Type>::offCentre_
(
    const GeoField& ddt0
) const
{
    // Multiply the ddt0 field by the spatially varying coefficient
    return ocCoeff_ * ddt0;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
template<class Type>
momentumConsistentCrankNicolson<Type>::momentumConsistentCrankNicolson
(
    const fvMesh& mesh,
    Istream& is
)
:
  ddtScheme<Type>(mesh,is),
  alpha_(mesh.lookupObject<volScalarField>("alpha.water")),
  bulkCnCoeff_(readScalar(is)),
  interfaceThreshold_(readScalar(is)),
  ocCoeff_
    (
        IOobject("ocCoeff", mesh.time().timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
        mesh,
        dimensionedScalar(dimless, bulkCnCoeff_)
    ),
  rDtCoefField_
    (
        IOobject("rDtCoef", mesh.time().timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
        mesh,
        dimensionedScalar(dimless/dimTime, 0.0)
    ),
  rDtCoef0Field_
    (
        IOobject("rDtCoef0", mesh.time().timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
        mesh,
        dimensionedScalar(dimless/dimTime, 0.0)
    )
{
    if (mesh.moving())
    {
        NotImplemented;
    }
}


template<class Type>
void momentumConsistentCrankNicolson<Type>::updateOcCoeff()
{
    const volVectorField& gradAlpha = fvc::grad(alpha_);
    
    forAll(ocCoeff_, cellI)
    {
        if (mag(gradAlpha[cellI]) > interfaceThreshold_)
        {
            ocCoeff_[cellI] = 0.0; // Euler
        }
        else
        {
            ocCoeff_[cellI] = bulkCnCoeff_; // CN
        }
    }
 
    ocCoeff_.correctBoundaryConditions();
    
    // --- DIAGNOSTIC START ---
    const bool debug=false;
    if (debug)
    {
        // FIX: Use tmp<> handle to avoid ambiguity
        tmp<surfaceScalarField> tOcSurf = fvc::interpolate(ocCoeff_);
        const surfaceScalarField& ocSurf = tOcSurf();
        
        const surfaceScalarField& magSf = mesh().magSf();
        const fvPatchList& patches = mesh().boundary();

        // 1. Internal Faces
        scalar globalArea = sum(magSf.primitiveField());
        scalar globalWeightedSum = sum(ocSurf.primitiveField() * magSf.primitiveField());
        scalar minFaceVal = gMin(ocSurf.primitiveField());
        scalar maxFaceVal = gMax(ocSurf.primitiveField());

        // 2. Boundary Faces
        forAll(ocSurf.boundaryField(), patchI)
        {
            const scalarField& pMagSf = magSf.boundaryField()[patchI];
            const scalarField& pOc = ocSurf.boundaryField()[patchI];

            // Calculate local sums for this patch
            scalar pArea = sum(pMagSf);
            scalar pWeighted = sum(pOc * pMagSf);
            
            // Check bounds on this patch
            if (pOc.size() > 0)
            {
                minFaceVal = min(minFaceVal, min(pOc));
                maxFaceVal = max(maxFaceVal, max(pOc));
            }

            // Handle Area Weighting for ALL coupled patches (Processor, Cyclic, etc.)
            if (patches[patchI].coupled())
            {
                globalArea += 0.5 * pArea;
                globalWeightedSum += 0.5 * pWeighted;
            }
            else
            {
                globalArea += pArea;
                globalWeightedSum += pWeighted;
            }
        }

        // Global Reduction
        reduce(globalArea, sumOp<scalar>());
        reduce(globalWeightedSum, sumOp<scalar>());
        reduce(minFaceVal, minOp<scalar>());
        reduce(maxFaceVal, maxOp<scalar>());

        if (Pstream::master())
        {
            Info<< "MomentumCN Diagnostic:" << nl
                << "  Face Avg = " << globalWeightedSum / (globalArea + VSMALL) << nl
                << "  Face Min = " << minFaceVal << " (Should be " << bulkCnCoeff_ << ")" << nl
                << "  Face Max = " << maxFaceVal << endl;
        }
    }
    // --- DIAGNOSTIC END ---

}

// =========================================================================
// fvmDdt Implementation
// =========================================================================
template<class Type>
tmp<fvMatrix<Type>>
momentumConsistentCrankNicolson<Type>::fvmDdt
(
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    DDt0Field<GeometricField<Type, fvPatchField, volMesh>>& ddt0 =
        ddt0_<GeometricField<Type, fvPatchField, volMesh>>
        (
            "ddt0(" + rho.name() + ',' + vf.name() + ')',
            rho.dimensions()*vf.dimensions()
        );

    tmp<fvMatrix<Type>> tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            rho.dimensions()*vf.dimensions()*dimVol/dimTime
        )
    );
    fvMatrix<Type>& fvm = tfvm.ref();

    fvm.diag() = rDtCoef_(ddt0)*rho.primitiveField()*mesh().V();

    vf.oldTime().oldTime();
    rho.oldTime().oldTime();
    
   if (mesh().moving())
    {
      NotImplemented;
    }
   else
     {
        if (evaluate(ddt0))
        {
            ddt0 = rDtCoef0_(ddt0)*
            (
                rho.oldTime()*vf.oldTime()
              - rho.oldTime().oldTime()*vf.oldTime().oldTime()
            ) - offCentre_(ddt0());
        }

        fvm.source() =
        (
	 rDtCoef_(ddt0)().primitiveField()*rho.oldTime().primitiveField()*vf.oldTime().primitiveField()
          + offCentre_(ddt0.primitiveField())
        )*mesh().V();
     }

   return tfvm;
}
  

// =========================================================================
// fvcDdtPhiCorr Implementation
// =========================================================================

template<class Type>
tmp<typename momentumConsistentCrankNicolson<Type>::fluxFieldType>
momentumConsistentCrankNicolson<Type>::fvcDdtPhiCorr
(
    const GeometricField<Type, fvPatchField, volMesh>& U,
    const fluxFieldType& phi
)
{
  //#include "gemini.H"
  #include "cn2.H"
  //NotImplemented; return nullptr;
}

// ... NotImplemented functions ...
template<class Type> tmp<GeometricField<Type, fvPatchField, volMesh>> momentumConsistentCrankNicolson<Type>::fvcDdt(const dimensioned<Type>& dt) { NotImplemented; return nullptr; }
template<class Type> tmp<GeometricField<Type, fvPatchField, volMesh>> momentumConsistentCrankNicolson<Type>::fvcDdt(const GeometricField<Type, fvPatchField, volMesh>& vf) { NotImplemented; return nullptr; }
template<class Type> tmp<GeometricField<Type, fvPatchField, volMesh>> momentumConsistentCrankNicolson<Type>::fvcDdt(const dimensionedScalar& rho, const GeometricField<Type, fvPatchField, volMesh>& vf) { NotImplemented; return nullptr; }
template<class Type> tmp<GeometricField<Type, fvPatchField, volMesh>> momentumConsistentCrankNicolson<Type>::fvcDdt(const volScalarField& rho, const GeometricField<Type, fvPatchField, volMesh>& vf) { NotImplemented; return nullptr; }
template<class Type> tmp<GeometricField<Type, fvPatchField, volMesh>> momentumConsistentCrankNicolson<Type>::fvcDdt(const volScalarField& alpha, const volScalarField& rho, const GeometricField<Type, fvPatchField, volMesh>& vf) { NotImplemented; return nullptr; }
template<class Type> tmp<fvMatrix<Type>> momentumConsistentCrankNicolson<Type>::fvmDdt(const GeometricField<Type, fvPatchField, volMesh>& vf) { NotImplemented; return nullptr; }
template<class Type> tmp<fvMatrix<Type>> momentumConsistentCrankNicolson<Type>::fvmDdt(const dimensionedScalar& rho, const GeometricField<Type, fvPatchField, volMesh>& vf) { NotImplemented; return nullptr; }
template<class Type> tmp<fvMatrix<Type>> momentumConsistentCrankNicolson<Type>::fvmDdt(const volScalarField& alpha, const volScalarField& rho, const GeometricField<Type, fvPatchField, volMesh>& vf) { NotImplemented; return nullptr; }
template<class Type> tmp<typename momentumConsistentCrankNicolson<Type>::fluxFieldType> momentumConsistentCrankNicolson<Type>::fvcDdtUfCorr(const GeometricField<Type, fvPatchField, volMesh>& U, const GeometricField<Type, fvsPatchField, surfaceMesh>& Uf) { NotImplemented; return fluxFieldType::null(); }
template<class Type> tmp<typename momentumConsistentCrankNicolson<Type>::fluxFieldType> momentumConsistentCrankNicolson<Type>::fvcDdtUfCorr(const volScalarField& rho, const GeometricField<Type, fvPatchField, volMesh>& U, const GeometricField<Type, fvsPatchField, surfaceMesh>& Uf) { NotImplemented; return fluxFieldType::null(); }
template<class Type> tmp<typename momentumConsistentCrankNicolson<Type>::fluxFieldType> momentumConsistentCrankNicolson<Type>::fvcDdtPhiCorr(const volScalarField& rho, const GeometricField<Type, fvPatchField, volMesh>& U, const fluxFieldType& phi) { NotImplemented; return fluxFieldType::null(); }
template<class Type> tmp<surfaceScalarField> momentumConsistentCrankNicolson<Type>::meshPhi(const GeometricField<Type, fvPatchField, volMesh>& vf) { NotImplemented; return nullptr; }

} // End namespace fv
} // End namespace Foam
