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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

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
    scalarField coeffs(mesh().nCells());

    if (mesh().time().timeIndex() > ddt0.startTimeIndex())
    {
        forAll(coeffs, cellI)
        {
            coeffs[cellI] = 1.0 + ocCoeff_[cellI];
        }
    }
    else
    {
        coeffs = 1.0;
    }

    return coeffs;
}


template<class Type>
template<class GeoField>
scalarField momentumConsistentCrankNicolson<Type>::coef0_
(
    const DDt0Field<GeoField>& ddt0
) const
{
    scalarField coeffs(mesh().nCells());

    if (mesh().time().timeIndex() > ddt0.startTimeIndex() + 1)
    {
        forAll(coeffs, cellI)
        {
            coeffs[cellI] = 1.0 + ocCoeff_[cellI];
        }
    }
    else
    {
        coeffs = 1.0;
    }

    return coeffs;
}


template<class Type>
template<class GeoField>
tmp<scalarField> momentumConsistentCrankNicolson<Type>::rDtCoef_
(
    const DDt0Field<GeoField>& ddt0
) const
{
    return coef_(ddt0)/mesh().time().deltaT().value();
}


template<class Type>
template<class GeoField>
tmp<scalarField> momentumConsistentCrankNicolson<Type>::rDtCoef0_
(
    const DDt0Field<GeoField>& ddt0
) const
{
    return coef0_(ddt0)/mesh().time().deltaT0().value();
}


template<class Type>
template<class GeoField>
tmp<GeoField> momentumConsistentCrankNicolson<Type>::offCentre_
(
    const GeoField& ddt0
) const
{
    tmp<GeoField> toffCentre
    (
        new GeoField
        (
            IOobject
            (
                "offCentre(" + ddt0.name() + ')',
                mesh().time().timeName(),
                mesh()
            ),
            mesh(),
            ddt0.dimensions()
        )
    );

    GeoField& offCentre = toffCentre.ref();

    forAll(offCentre.primitiveField(), cellI)
    {
        offCentre.primitiveFieldRef()[cellI] =
            ocCoeff_[cellI]*ddt0.primitiveField()[cellI];
    }

    offCentre.boundaryFieldRef() = ddt0.boundaryField();

    return toffCentre;
}


// * * * * * * * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
momentumConsistentCrankNicolson<Type>::momentumConsistentCrankNicolson
(
    const fvMesh& mesh,
    Istream& is
)
:
    ddtScheme<Type>(mesh),
    alpha_(mesh.lookupObject<volScalarField>("alpha.water")),
    bulkCnCoeff_(readScalar(is)),
    interfaceThreshold_(readScalar(is)),
    ocCoeff_(mesh.nCells(), bulkCnCoeff_)
{
    if (mesh.moving())
    {
        mesh.V00();
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void momentumConsistentCrankNicolson<Type>::updateOcCoeff() const
{

    tmp<volVectorField> tgradAlpha = fvc::grad(alpha_);
    const volVectorField& gradAlpha = tgradAlpha();

    forAll(ocCoeff_, cellI)
    {
        scalar gradMag = mag(gradAlpha[cellI]);

        if (gradMag > interfaceThreshold_)
        {
            ocCoeff_[cellI] = 0.0;
        }
        else
        {
            ocCoeff_[cellI] = bulkCnCoeff_;
        }
    }
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
momentumConsistentCrankNicolson<Type>::fvcDdt
(
    const dimensioned<Type>& dt
)
{
    DDt0Field<GeometricField<Type, fvPatchField, volMesh>>& ddt0 =
        ddt0_<GeometricField<Type, fvPatchField, volMesh>>
        (
            "ddt0(" + dt.name() + ')',
            dt.dimensions()
        );

    IOobject ddtIOobject
    (
        "ddt(" + dt.name() + ')',
        mesh().time().timeName(),
        mesh().thisDb()
    );

    tmp<GeometricField<Type, fvPatchField, volMesh>> tdtdt
    (
        new GeometricField<Type, fvPatchField, volMesh>
        (
            ddtIOobject,
            mesh(),
            Foam::zero{},
            (dt.dimensions()/dimTime)
        )
    );

    updateOcCoeff();
    tmp<scalarField> rDtCoef = rDtCoef_(ddt0);

    if (mesh().moving())
    {
        if (evaluate(ddt0))
        {
            tmp<scalarField> rDtCoef0 = rDtCoef0_(ddt0);

            forAll(ddt0.primitiveFieldRef(), cellI)
            {
                ddt0.primitiveFieldRef()[cellI] =
                (
                    (rDtCoef0()[cellI]*dt.value())
                   *(mesh().V0()[cellI] - mesh().V00()[cellI])
                  - mesh().V00()[cellI]
                   *ocCoeff_[cellI]*ddt0.primitiveField()[cellI]
                )/mesh().V0()[cellI];
            }
        }

        forAll(tdtdt.ref().primitiveFieldRef(), cellI)
        {
            tdtdt.ref().primitiveFieldRef()[cellI] =
            (
                (rDtCoef()[cellI]*dt.value())
               *(mesh().V()[cellI] - mesh().V0()[cellI])
              - mesh().V0()[cellI]
               *ocCoeff_[cellI]*ddt0.primitiveField()[cellI]
            )/mesh().V()[cellI];
        }

        tdtdt.ref().boundaryFieldRef().
            template evaluateCoupled<coupledFvPatch>();
    }

    return tdtdt;
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
momentumConsistentCrankNicolson<Type>::fvcDdt
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    DDt0Field<GeometricField<Type, fvPatchField, volMesh>>& ddt0 =
        ddt0_<GeometricField<Type, fvPatchField, volMesh>>
        (
            "ddt0(" + vf.name() + ')',
            vf.dimensions()
        );

    IOobject ddtIOobject
    (
        "ddt(" + vf.name() + ')',
        mesh().time().timeName(),
        mesh().thisDb()
    );

    updateOcCoeff();
    tmp<scalarField> trDtCoef = rDtCoef_(ddt0);
    const scalarField& rDtCoef = trDtCoef();

    if (mesh().moving())
    {
        if (evaluate(ddt0))
        {
            tmp<scalarField> trDtCoef0 = rDtCoef0_(ddt0);
            const scalarField& rDtCoef0 = trDtCoef0();

            forAll(ddt0.primitiveFieldRef(), cellI)
            {
                ddt0.primitiveFieldRef()[cellI] =
                (
                    rDtCoef0[cellI]*
                    (
                        mesh().V0()[cellI]*vf.oldTime().primitiveField()[cellI]
                      - mesh().V00()[cellI]
                       *vf.oldTime().oldTime().primitiveField()[cellI]
                    )
                  - mesh().V00()[cellI]
                   *ocCoeff_[cellI]*ddt0.primitiveField()[cellI]
                )/mesh().V0()[cellI];
            }

            typename GeometricField<Type, fvPatchField, volMesh>::Boundary& ddt0bf = ddt0.boundaryFieldRef();
            const GeometricField<Type, fvPatchField, volMesh>& vfOld = vf.oldTime();
            const GeometricField<Type, fvPatchField, volMesh>& vfOldOld = vf.oldTime().oldTime();
            forAll(ddt0bf, patchI)
            {
                fvPatchField<Type>& pddt0 = ddt0bf[patchI];
                const fvPatchField<Type>& pvfOld = vfOld.boundaryField()[patchI];
                const fvPatchField<Type>& pvfOldOld = vfOldOld.boundaryField()[patchI];
                const labelUList& faceCells = pddt0.patch().faceCells();
                forAll(pddt0, i)
                {
                    label own = faceCells[i];
                    pddt0[i] = rDtCoef0[own] * (pvfOld[i] - pvfOldOld[i]);
                }
            }
        }

        tmp<GeometricField<Type, fvPatchField, volMesh>> tdtdt
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                mesh(),
                vf.dimensions()/dimTime
            )
        );

        forAll(tdtdt.ref().primitiveFieldRef(), cellI)
        {
            tdtdt.ref().primitiveFieldRef()[cellI] =
            (
                rDtCoef[cellI]*
                (
                    mesh().V()[cellI]*vf.primitiveField()[cellI]
                  - mesh().V0()[cellI]*vf.oldTime().primitiveField()[cellI]
                )
              - mesh().V0()[cellI]
               *ocCoeff_[cellI]*ddt0.primitiveField()[cellI]
            )/mesh().V()[cellI];
        }

        typename GeometricField<Type, fvPatchField, volMesh>::Boundary& dtdtbf = tdtdt.ref().boundaryFieldRef();
        const GeometricField<Type, fvPatchField, volMesh>& vfCurrent = vf;
        const GeometricField<Type, fvPatchField, volMesh>& vfOld = vf.oldTime();
        forAll(dtdtbf, patchI)
        {
            fvPatchField<Type>& pdtdt = dtdtbf[patchI];
            const fvPatchField<Type>& pvfCurrent = vfCurrent.boundaryField()[patchI];
            const fvPatchField<Type>& pvfOld = vfOld.boundaryField()[patchI];
            const labelUList& faceCells = pdtdt.patch().faceCells();
            forAll(pdtdt, i)
            {
                label own = faceCells[i];
                pdtdt[i] = rDtCoef[own] * (pvfCurrent[i] - pvfOld[i]);
            }
        }

        tdtdt.ref().boundaryFieldRef().
            template evaluateCoupled<coupledFvPatch>();

        return tdtdt;
    }
    else
    {
        if (evaluate(ddt0))
        {
            tmp<scalarField> trDtCoef0 = rDtCoef0_(ddt0);
            const scalarField& rDtCoef0 = trDtCoef0();

            forAll(ddt0.primitiveFieldRef(), cellI)
            {
                ddt0.primitiveFieldRef()[cellI] =
                    rDtCoef0[cellI]*
                    (
                        vf.oldTime().primitiveField()[cellI]
                      - vf.oldTime().oldTime().primitiveField()[cellI]
                    )
                  - ocCoeff_[cellI]*ddt0.primitiveField()[cellI];
            }

            ddt0.boundaryFieldRef() = vf.oldTime().boundaryField();
        }

        tmp<GeometricField<Type, fvPatchField, volMesh>> tdtdt
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                mesh(),
                vf.dimensions()/dimTime
            )
        );

        forAll(tdtdt.ref().primitiveFieldRef(), cellI)
        {
            tdtdt.ref().primitiveFieldRef()[cellI] =
                rDtCoef[cellI]*
                (
                    vf.primitiveField()[cellI]
                  - vf.oldTime().primitiveField()[cellI]
                )
              - ocCoeff_[cellI]*ddt0.primitiveField()[cellI];
        }

        tdtdt.ref().boundaryFieldRef() = vf.boundaryField();

        return tdtdt;
    }
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
momentumConsistentCrankNicolson<Type>::fvcDdt
(
    const dimensionedScalar& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    DDt0Field<GeometricField<Type, fvPatchField, volMesh>>& ddt0 =
        ddt0_<GeometricField<Type, fvPatchField, volMesh>>
        (
            "ddt0(" + rho.name() + ',' + vf.name() + ')',
            rho.dimensions()*vf.dimensions()
        );

    IOobject ddtIOobject
    (
        "ddt(" + rho.name() + ',' + vf.name() + ')',
        mesh().time().timeName(),
        mesh().thisDb()
    );

    updateOcCoeff();
    tmp<scalarField> trDtCoef = rDtCoef_(ddt0);
    const scalarField& rDtCoef = trDtCoef();

    if (mesh().moving())
    {
        if (evaluate(ddt0))
        {
            tmp<scalarField> trDtCoef0 = rDtCoef0_(ddt0);
            const scalarField& rDtCoef0 = trDtCoef0();

            forAll(ddt0.primitiveFieldRef(), cellI)
            {
                ddt0.primitiveFieldRef()[cellI] =
                (
                    rDtCoef0[cellI]*rho.value()*
                    (
                        mesh().V0()[cellI]*vf.oldTime().primitiveField()[cellI]
                      - mesh().V00()[cellI]
                       *vf.oldTime().oldTime().primitiveField()[cellI]
                    )
                  - mesh().V00()[cellI]
                   *ocCoeff_[cellI]*ddt0.primitiveField()[cellI]
                )/mesh().V0()[cellI];
            }

            typename GeometricField<Type, fvPatchField, volMesh>::Boundary& ddt0bf = ddt0.boundaryFieldRef();
            const GeometricField<Type, fvPatchField, volMesh>& vfOld = vf.oldTime();
            const GeometricField<Type, fvPatchField, volMesh>& vfOldOld = vf.oldTime().oldTime();
            forAll(ddt0bf, patchI)
            {
                fvPatchField<Type>& pddt0 = ddt0bf[patchI];
                const fvPatchField<Type>& pvfOld = vfOld.boundaryField()[patchI];
                const fvPatchField<Type>& pvfOldOld = vfOldOld.boundaryField()[patchI];
                const labelUList& faceCells = pddt0.patch().faceCells();
                forAll(pddt0, i)
                {
                    label own = faceCells[i];
                    pddt0[i] = rDtCoef0[own] * rho.value() * (pvfOld[i] - pvfOldOld[i]);
                }
            }
        }

        tmp<GeometricField<Type, fvPatchField, volMesh>> tdtdt
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                mesh(),
                rho.dimensions()*vf.dimensions()/dimTime
            )
        );

        forAll(tdtdt.ref().primitiveFieldRef(), cellI)
        {
            tdtdt.ref().primitiveFieldRef()[cellI] =
            (
                rDtCoef[cellI]*rho.value()*
                (
                    mesh().V()[cellI]*vf.primitiveField()[cellI]
                  - mesh().V0()[cellI]*vf.oldTime().primitiveField()[cellI]
                )
              - mesh().V0()[cellI]
               *ocCoeff_[cellI]*ddt0.primitiveField()[cellI]
            )/mesh().V()[cellI];
        }

        typename GeometricField<Type, fvPatchField, volMesh>::Boundary& dtdtbf = tdtdt.ref().boundaryFieldRef();
        const GeometricField<Type, fvPatchField, volMesh>& vfCurrent = vf;
        const GeometricField<Type, fvPatchField, volMesh>& vfOld = vf.oldTime();
        forAll(dtdtbf, patchI)
        {
            fvPatchField<Type>& pdtdt = dtdtbf[patchI];
            const fvPatchField<Type>& pvfCurrent = vfCurrent.boundaryField()[patchI];
            const fvPatchField<Type>& pvfOld = vfOld.boundaryField()[patchI];
            const labelUList& faceCells = pdtdt.patch().faceCells();
            forAll(pdtdt, i)
            {
                label own = faceCells[i];
                pdtdt[i] = rDtCoef[own] * rho.value() * (pvfCurrent[i] - pvfOld[i]);
            }
        }

        tdtdt.ref().boundaryFieldRef().
            template evaluateCoupled<coupledFvPatch>();

        return tdtdt;
    }
    else
    {
        if (evaluate(ddt0))
        {
            tmp<scalarField> trDtCoef0 = rDtCoef0_(ddt0);
            const scalarField& rDtCoef0 = trDtCoef0();

            forAll(ddt0.primitiveFieldRef(), cellI)
            {
                ddt0.primitiveFieldRef()[cellI] =
                    rDtCoef0[cellI]*rho.value()*
                    (
                        vf.oldTime().primitiveField()[cellI]
                      - vf.oldTime().oldTime().primitiveField()[cellI]
                    )
                  - ocCoeff_[cellI]*ddt0.primitiveField()[cellI];
            }
        }

        tmp<GeometricField<Type, fvPatchField, volMesh>> tdtdt
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                mesh(),
                rho.dimensions()*vf.dimensions()/dimTime
            )
        );

        forAll(tdtdt.ref().primitiveFieldRef(), cellI)
        {
            tdtdt.ref().primitiveFieldRef()[cellI] =
                rDtCoef[cellI]*rho.value()*
                (
                    vf.primitiveField()[cellI]
                  - vf.oldTime().primitiveField()[cellI]
                )
              - ocCoeff_[cellI]*ddt0.primitiveField()[cellI];
        }

        return tdtdt;
    }
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
momentumConsistentCrankNicolson<Type>::fvcDdt
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

    IOobject ddtIOobject
    (
        "ddt(" + rho.name() + ',' + vf.name() + ')',
        mesh().time().timeName(),
        mesh().thisDb()
    );

    updateOcCoeff();
    tmp<scalarField> trDtCoef = rDtCoef_(ddt0);
    const scalarField& rDtCoef = trDtCoef();

    if (mesh().moving())
    {
        if (evaluate(ddt0))
        {
            tmp<scalarField> trDtCoef0 = rDtCoef0_(ddt0);
            const scalarField& rDtCoef0 = trDtCoef0();

            forAll(ddt0.primitiveFieldRef(), cellI)
            {
                ddt0.primitiveFieldRef()[cellI] =
                (
                    rDtCoef0[cellI]*
                    (
                        mesh().V0()[cellI]
                       *rho.oldTime().primitiveField()[cellI]
                       *vf.oldTime().primitiveField()[cellI]
                      - mesh().V00()[cellI]
                       *rho.oldTime().oldTime().primitiveField()[cellI]
                       *vf.oldTime().oldTime().primitiveField()[cellI]
                    )
                  - mesh().V00()[cellI]
                   *ocCoeff_[cellI]*ddt0.primitiveField()[cellI]
                )/mesh().V0()[cellI];
            }

            typename GeometricField<Type, fvPatchField, volMesh>::Boundary& ddt0bf = ddt0.boundaryFieldRef();
            const volScalarField& rhoOld = rho.oldTime();
            const volScalarField& rhoOldOld = rho.oldTime().oldTime();
            const GeometricField<Type, fvPatchField, volMesh>& vfOld = vf.oldTime();
            const GeometricField<Type, fvPatchField, volMesh>& vfOldOld = vf.oldTime().oldTime();
            forAll(ddt0bf, patchI)
            {
                fvPatchField<Type>& pddt0 = ddt0bf[patchI];
                const fvPatchScalarField& prhoOld = rhoOld.boundaryField()[patchI];
                const fvPatchScalarField& prhoOldOld = rhoOldOld.boundaryField()[patchI];
                const fvPatchField<Type>& pvfOld = vfOld.boundaryField()[patchI];
                const fvPatchField<Type>& pvfOldOld = vfOldOld.boundaryField()[patchI];
                const labelUList& faceCells = pddt0.patch().faceCells();
                forAll(pddt0, i)
                {
                    label own = faceCells[i];
                    pddt0[i] = rDtCoef0[own] * (prhoOld[i] * pvfOld[i] - prhoOldOld[i] * pvfOldOld[i]);
                }
            }
        }

        tmp<GeometricField<Type, fvPatchField, volMesh>> tdtdt
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                mesh(),
                rho.dimensions()*vf.dimensions()/dimTime
            )
        );

        forAll(tdtdt.ref().primitiveFieldRef(), cellI)
        {
            tdtdt.ref().primitiveFieldRef()[cellI] =
            (
                rDtCoef[cellI]*
                (
                    mesh().V()[cellI]
                   *rho.primitiveField()[cellI]
                   *vf.primitiveField()[cellI]
                  - mesh().V0()[cellI]
                   *rho.oldTime().primitiveField()[cellI]
                   *vf.oldTime().primitiveField()[cellI]
                )
              - mesh().V0()[cellI]
               *ocCoeff_[cellI]*ddt0.primitiveField()[cellI]
            )/mesh().V()[cellI];
        }

        typename GeometricField<Type, fvPatchField, volMesh>::Boundary& dtdtbf = tdtdt.ref().boundaryFieldRef();
        const volScalarField& rhoCurrent = rho;
        const volScalarField& rhoOld = rho.oldTime();
        const GeometricField<Type, fvPatchField, volMesh>& vfCurrent = vf;
        const GeometricField<Type, fvPatchField, volMesh>& vfOld = vf.oldTime();
        forAll(dtdtbf, patchI)
        {
            fvPatchField<Type>& pdtdt = dtdtbf[patchI];
            const fvPatchScalarField& prhoCurrent = rhoCurrent.boundaryField()[patchI];
            const fvPatchScalarField& prhoOld = rhoOld.boundaryField()[patchI];
            const fvPatchField<Type>& pvfCurrent = vfCurrent.boundaryField()[patchI];
            const fvPatchField<Type>& pvfOld = vfOld.boundaryField()[patchI];
            const labelUList& faceCells = pdtdt.patch().faceCells();
            forAll(pdtdt, i)
            {
                label own = faceCells[i];
                pdtdt[i] = rDtCoef[own] * (prhoCurrent[i] * pvfCurrent[i] - prhoOld[i] * pvfOld[i]);
            }
        }

        tdtdt.ref().boundaryFieldRef().
            template evaluateCoupled<coupledFvPatch>();

        return tdtdt;
    }
    else
    {
        if (evaluate(ddt0))
        {
            tmp<scalarField> trDtCoef0 = rDtCoef0_(ddt0);
            const scalarField& rDtCoef0 = trDtCoef0();

            forAll(ddt0.primitiveFieldRef(), cellI)
            {
                ddt0.primitiveFieldRef()[cellI] =
                    rDtCoef0[cellI]*
                    (
                        rho.oldTime().primitiveField()[cellI]
                       *vf.oldTime().primitiveField()[cellI]
                      - rho.oldTime().oldTime().primitiveField()[cellI]
                       *vf.oldTime().oldTime().primitiveField()[cellI]
                    )
                  - ocCoeff_[cellI]*ddt0.primitiveField()[cellI];
            }
        }

        tmp<GeometricField<Type, fvPatchField, volMesh>> tdtdt
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                mesh(),
                rho.dimensions()*vf.dimensions()/dimTime
            )
        );

        forAll(tdtdt.ref().primitiveFieldRef(), cellI)
        {
            tdtdt.ref().primitiveFieldRef()[cellI] =
                rDtCoef[cellI]*
                (
                    rho.primitiveField()[cellI]*vf.primitiveField()[cellI]
                  - rho.oldTime().primitiveField()[cellI]
                   *vf.oldTime().primitiveField()[cellI]
                )
              - ocCoeff_[cellI]*ddt0.primitiveField()[cellI];
        }

        return tdtdt;
    }
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
momentumConsistentCrankNicolson<Type>::fvcDdt
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    DDt0Field<GeometricField<Type, fvPatchField, volMesh>>& ddt0 =
        ddt0_<GeometricField<Type, fvPatchField, volMesh>>
        (
            "ddt0(" + alpha.name() + ',' + rho.name() + ',' + vf.name() + ')',
            alpha.dimensions()*rho.dimensions()*vf.dimensions()
        );

    IOobject ddtIOobject
    (
        "ddt(" + alpha.name() + ',' + rho.name() + ',' + vf.name() + ')',
        mesh().time().timeName(),
        mesh().thisDb()
    );

    updateOcCoeff();
    tmp<scalarField> trDtCoef = rDtCoef_(ddt0);
    const scalarField& rDtCoef = trDtCoef();

    if (mesh().moving())
    {
        if (evaluate(ddt0))
        {
            tmp<scalarField> trDtCoef0 = rDtCoef0_(ddt0);
            const scalarField& rDtCoef0 = trDtCoef0();

            forAll(ddt0.primitiveFieldRef(), cellI)
            {
                ddt0.primitiveFieldRef()[cellI] =
                (
                    rDtCoef0[cellI]*
                    (
                        mesh().V0()[cellI]
                       *alpha.oldTime().primitiveField()[cellI]
                       *rho.oldTime().primitiveField()[cellI]
                       *vf.oldTime().primitiveField()[cellI]

                      - mesh().V00()[cellI]
                       *alpha.oldTime().oldTime().primitiveField()[cellI]
                       *rho.oldTime().oldTime().primitiveField()[cellI]
                       *vf.oldTime().oldTime().primitiveField()[cellI]
                    )
                  - mesh().V00()[cellI]
                   *ocCoeff_[cellI]*ddt0.primitiveField()[cellI]
                )/mesh().V0()[cellI];
            }

            typename GeometricField<Type, fvPatchField, volMesh>::Boundary& ddt0bf = ddt0.boundaryFieldRef();
            const volScalarField& alphaOld = alpha.oldTime();
            const volScalarField& alphaOldOld = alpha.oldTime().oldTime();
            const volScalarField& rhoOld = rho.oldTime();
            const volScalarField& rhoOldOld = rho.oldTime().oldTime();
            const GeometricField<Type, fvPatchField, volMesh>& vfOld = vf.oldTime();
            const GeometricField<Type, fvPatchField, volMesh>& vfOldOld = vf.oldTime().oldTime();
            forAll(ddt0bf, patchI)
            {
                fvPatchField<Type>& pddt0 = ddt0bf[patchI];
                const fvPatchScalarField& palphaOld = alphaOld.boundaryField()[patchI];
                const fvPatchScalarField& palphaOldOld = alphaOldOld.boundaryField()[patchI];
                const fvPatchScalarField& prhoOld = rhoOld.boundaryField()[patchI];
                const fvPatchScalarField& prhoOldOld = rhoOldOld.boundaryField()[patchI];
                const fvPatchField<Type>& pvfOld = vfOld.boundaryField()[patchI];
                const fvPatchField<Type>& pvfOldOld = vfOldOld.boundaryField()[patchI];
                const labelUList& faceCells = pddt0.patch().faceCells();
                forAll(pddt0, i)
                {
                    label own = faceCells[i];
                    pddt0[i] = rDtCoef0[own] *
                               (palphaOld[i] * prhoOld[i] * pvfOld[i] -
                                palphaOldOld[i] * prhoOldOld[i] * pvfOldOld[i]);
                }
            }
        }

        tmp<GeometricField<Type, fvPatchField, volMesh>> tdtdt
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                mesh(),
               alpha.dimensions()*rho.dimensions()*vf.dimensions()/dimTime
            )
        );

        forAll(tdtdt.ref().primitiveFieldRef(), cellI)
        {
            tdtdt.ref().primitiveFieldRef()[cellI] =
            (
                rDtCoef[cellI]*
                (
                    mesh().V()[cellI]
                   *alpha.primitiveField()[cellI]
                   *rho.primitiveField()[cellI]
                   *vf.primitiveField()[cellI]

                  - mesh().V0()[cellI]
                   *alpha.oldTime().primitiveField()[cellI]
                   *rho.oldTime().primitiveField()[cellI]
                   *vf.oldTime().primitiveField()[cellI]
                )
              - mesh().V0()[cellI]
               *ocCoeff_[cellI]*ddt0.primitiveField()[cellI]
            )/mesh().V()[cellI];
        }

        typename GeometricField<Type, fvPatchField, volMesh>::Boundary& dtdtbf = tdtdt.ref().boundaryFieldRef();
        const volScalarField& alphaCurrent = alpha;
        const volScalarField& alphaOld = alpha.oldTime();
        const volScalarField& rhoCurrent = rho;
        const volScalarField& rhoOld = rho.oldTime();
        const GeometricField<Type, fvPatchField, volMesh>& vfCurrent = vf;
        const GeometricField<Type, fvPatchField, volMesh>& vfOld = vf.oldTime();
        forAll(dtdtbf, patchI)
        {
            fvPatchField<Type>& pdtdt = dtdtbf[patchI];
            const fvPatchScalarField& palphaCurrent = alphaCurrent.boundaryField()[patchI];
            const fvPatchScalarField& palphaOld = alphaOld.boundaryField()[patchI];
            const fvPatchScalarField& prhoCurrent = rhoCurrent.boundaryField()[patchI];
            const fvPatchScalarField& prhoOld = rhoOld.boundaryField()[patchI];
            const fvPatchField<Type>& pvfCurrent = vfCurrent.boundaryField()[patchI];
            const fvPatchField<Type>& pvfOld = vfOld.boundaryField()[patchI];
            const labelUList& faceCells = pdtdt.patch().faceCells();
            forAll(pdtdt, i)
            {
                label own = faceCells[i];
                pdtdt[i] = rDtCoef[own] *
                           (palphaCurrent[i] * prhoCurrent[i] * pvfCurrent[i] -
                            palphaOld[i] * prhoOld[i] * pvfOld[i]);
            }
        }

        tdtdt.ref().boundaryFieldRef().
            template evaluateCoupled<coupledFvPatch>();

        return tdtdt;
    }
    else
    {
        if (evaluate(ddt0))
        {
            tmp<scalarField> trDtCoef0 = rDtCoef0_(ddt0);
            const scalarField& rDtCoef0 = trDtCoef0();

            forAll(ddt0.primitiveFieldRef(), cellI)
            {
                ddt0.primitiveFieldRef()[cellI] =
                    rDtCoef0[cellI]*
                    (
                        alpha.oldTime().primitiveField()[cellI]
                       *rho.oldTime().primitiveField()[cellI]
                       *vf.oldTime().primitiveField()[cellI]

                      - alpha.oldTime().oldTime().primitiveField()[cellI]
                       *rho.oldTime().oldTime().primitiveField()[cellI]
                       *vf.oldTime().oldTime().primitiveField()[cellI]
                    )
                  - ocCoeff_[cellI]*ddt0.primitiveField()[cellI];
            }
        }

        tmp<GeometricField<Type, fvPatchField, volMesh>> tdtdt
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                mesh(),
               alpha.dimensions()*rho.dimensions()*vf.dimensions()/dimTime
            )
        );

        forAll(tdtdt.ref().primitiveFieldRef(), cellI)
        {
            tdtdt.ref().primitiveFieldRef()[cellI] =
                rDtCoef[cellI]*
                (
                    alpha.primitiveField()[cellI]
                   *rho.primitiveField()[cellI]
                   *vf.primitiveField()[cellI]
                  - alpha.oldTime().primitiveField()[cellI]
                   *rho.oldTime().primitiveField()[cellI]
                   *vf.oldTime().primitiveField()[cellI]
                )
              - ocCoeff_[cellI]*ddt0.primitiveField()[cellI];
        }

        return tdtdt;
    }
}


template<class Type>
tmp<fvMatrix<Type>>
momentumConsistentCrankNicolson<Type>::fvmDdt
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    DDt0Field<GeometricField<Type, fvPatchField, volMesh>>& ddt0 =
        ddt0_<GeometricField<Type, fvPatchField, volMesh>>
        (
            "ddt0(" + vf.name() + ')',
            vf.dimensions()
        );

    tmp<fvMatrix<Type>> tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            vf.dimensions()*dimVol/dimTime
        )
    );

    fvMatrix<Type>& fvm = tfvm.ref();

    updateOcCoeff();
    tmp<scalarField> trDtCoef = rDtCoef_(ddt0);
    const scalarField& rDtCoef = trDtCoef();

    const scalarField& V = mesh().V();
    scalarField& diag = fvm.diag();
    Field<Type>& source = fvm.source();

    forAll(diag, cellI)
    {
        diag[cellI] = rDtCoef[cellI]*V[cellI];
    }

    vf.oldTime().oldTime();

    if (mesh().moving())
    {
        if (evaluate(ddt0))
        {
            tmp<scalarField> trDtCoef0 = rDtCoef0_(ddt0);
            const scalarField& rDtCoef0 = trDtCoef0();

            forAll(ddt0.primitiveFieldRef(), cellI)
            {
                ddt0.primitiveFieldRef()[cellI] =
                (
                    rDtCoef0[cellI]*
                    (
                        mesh().V0()[cellI]*vf.oldTime().primitiveField()[cellI]
                      - mesh().V00()[cellI]
                       *vf.oldTime().oldTime().primitiveField()[cellI]
                    )
                  - mesh().V00()[cellI]
                   *ocCoeff_[cellI]*ddt0.primitiveField()[cellI]
                )/mesh().V0()[cellI];
            }

            typename GeometricField<Type, fvPatchField, volMesh>::Boundary& ddt0bf = ddt0.boundaryFieldRef();
            const GeometricField<Type, fvPatchField, volMesh>& vfOld = vf.oldTime();
            const GeometricField<Type, fvPatchField, volMesh>& vfOldOld = vf.oldTime().oldTime();
            forAll(ddt0bf, patchI)
            {
                fvPatchField<Type>& pddt0 = ddt0bf[patchI];
                const fvPatchField<Type>& pvfOld = vfOld.boundaryField()[patchI];
                const fvPatchField<Type>& pvfOldOld = vfOldOld.boundaryField()[patchI];
                const labelUList& faceCells = pddt0.patch().faceCells();
                forAll(pddt0, i)
                {
                    label own = faceCells[i];
                    pddt0[i] = rDtCoef0[own] * (pvfOld[i] - pvfOldOld[i]);
                }
            }
        }

        forAll(source, cellI)
        {
            source[cellI] =
            (
                rDtCoef[cellI]*vf.oldTime().primitiveField()[cellI]
              + ocCoeff_[cellI]*ddt0.primitiveField()[cellI]
            )*mesh().V0()[cellI];
        }
    }
    else
    {
        if (evaluate(ddt0))
        {
            tmp<scalarField> trDtCoef0 = rDtCoef0_(ddt0);
            const scalarField& rDtCoef0 = trDtCoef0();

            forAll(ddt0.primitiveFieldRef(), cellI)
            {
                ddt0.primitiveFieldRef()[cellI] =
                    rDtCoef0[cellI]*
                    (
                        vf.oldTime().primitiveField()[cellI]
                      - vf.oldTime().oldTime().primitiveField()[cellI]
                    )
                  - ocCoeff_[cellI]*ddt0.primitiveField()[cellI];
            }
        }

        forAll(source, cellI)
        {
            source[cellI] =
            (
                rDtCoef[cellI]*vf.oldTime().primitiveField()[cellI]
              + ocCoeff_[cellI]*ddt0.primitiveField()[cellI]
            )*V[cellI];
        }
    }

    return tfvm;
}


template<class Type>
tmp<fvMatrix<Type>>
momentumConsistentCrankNicolson<Type>::fvmDdt
(
    const dimensionedScalar& rho,
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

    updateOcCoeff();
    tmp<scalarField> trDtCoef = rDtCoef_(ddt0);
    const scalarField& rDtCoef = trDtCoef();

    const scalarField& V = mesh().V();
    scalarField& diag = fvm.diag();
    Field<Type>& source = fvm.source();

    forAll(diag, cellI)
    {
        diag[cellI] = rDtCoef[cellI]*rho.value()*V[cellI];
    }

    vf.oldTime().oldTime();

    if (mesh().moving())
    {
        if (evaluate(ddt0))
        {
            tmp<scalarField> trDtCoef0 = rDtCoef0_(ddt0);
            const scalarField& rDtCoef0 = trDtCoef0();

            forAll(ddt0.primitiveFieldRef(), cellI)
            {
                ddt0.primitiveFieldRef()[cellI] =
                (
                    rDtCoef0[cellI]*rho.value()*
                    (
                        mesh().V0()[cellI]*vf.oldTime().primitiveField()[cellI]
                      - mesh().V00()[cellI]
                       *vf.oldTime().oldTime().primitiveField()[cellI]
                    )
                  - mesh().V00()[cellI]
                   *ocCoeff_[cellI]*ddt0.primitiveField()[cellI]
                )/mesh().V0()[cellI];
            }

            typename GeometricField<Type, fvPatchField, volMesh>::Boundary& ddt0bf = ddt0.boundaryFieldRef();
            const GeometricField<Type, fvPatchField, volMesh>& vfOld = vf.oldTime();
            const GeometricField<Type, fvPatchField, volMesh>& vfOldOld = vf.oldTime().oldTime();
            forAll(ddt0bf, patchI)
            {
                fvPatchField<Type>& pddt0 = ddt0bf[patchI];
                const fvPatchField<Type>& pvfOld = vfOld.boundaryField()[patchI];
                const fvPatchField<Type>& pvfOldOld = vfOldOld.boundaryField()[patchI];
                const labelUList& faceCells = pddt0.patch().faceCells();
                forAll(pddt0, i)
                {
                    label own = faceCells[i];
                    pddt0[i] = rDtCoef0[own] * rho.value() * (pvfOld[i] - pvfOldOld[i]);
                }
            }
        }

        forAll(source, cellI)
        {
            source[cellI] =
            (
                rDtCoef[cellI]*rho.value()*vf.oldTime().primitiveField()[cellI]
              + ocCoeff_[cellI]*ddt0.primitiveField()[cellI]
            )*mesh().V0()[cellI];
        }
    }
    else
    {
        if (evaluate(ddt0))
        {
            tmp<scalarField> trDtCoef0 = rDtCoef0_(ddt0);
            const scalarField& rDtCoef0 = trDtCoef0();

            forAll(ddt0.primitiveFieldRef(), cellI)
            {
                ddt0.primitiveFieldRef()[cellI] =
                    rDtCoef0[cellI]*rho.value()*
                    (
                        vf.oldTime().primitiveField()[cellI]
                      - vf.oldTime().oldTime().primitiveField()[cellI]
                    )
                  - ocCoeff_[cellI]*ddt0.primitiveField()[cellI];
            }
        }

        forAll(source, cellI)
        {
            source[cellI] =
            (
                rDtCoef[cellI]*rho.value()*vf.oldTime().primitiveField()[cellI]
              + ocCoeff_[cellI]*ddt0.primitiveField()[cellI]
            )*V[cellI];
        }
    }

    return tfvm;
}


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

    updateOcCoeff();
    tmp<scalarField> trDtCoef = rDtCoef_(ddt0);
    const scalarField& rDtCoef = trDtCoef();

    const scalarField& V = mesh().V();
    const scalarField& rhoField = rho.primitiveField();
    scalarField& diag = fvm.diag();
    Field<Type>& source = fvm.source();

    forAll(diag, cellI)
    {
        diag[cellI] = rDtCoef[cellI]*rhoField[cellI]*V[cellI];
    }

    vf.oldTime().oldTime();
    rho.oldTime().oldTime();

    if (mesh().moving())
    {
        if (evaluate(ddt0))
        {
            tmp<scalarField> trDtCoef0 = rDtCoef0_(ddt0);
            const scalarField& rDtCoef0 = trDtCoef0();

            forAll(ddt0.primitiveFieldRef(), cellI)
            {
                ddt0.primitiveFieldRef()[cellI] =
                (
                    rDtCoef0[cellI]*
                    (
                        mesh().V0()[cellI]*rho.oldTime().primitiveField()[cellI]
                       *vf.oldTime().primitiveField()[cellI]
                      - mesh().V00()[cellI]*rho.oldTime().oldTime().primitiveField()[cellI]
                       *vf.oldTime().oldTime().primitiveField()[cellI]
                    )
                  - mesh().V00()[cellI]
                   *ocCoeff_[cellI]*ddt0.primitiveField()[cellI]
                )/mesh().V0()[cellI];
            }

            typename GeometricField<Type, fvPatchField, volMesh>::Boundary& ddt0bf = ddt0.boundaryFieldRef();
            const volScalarField& rhoOld = rho.oldTime();
            const volScalarField& rhoOldOld = rho.oldTime().oldTime();
            const GeometricField<Type, fvPatchField, volMesh>& vfOld = vf.oldTime();
            const GeometricField<Type, fvPatchField, volMesh>& vfOldOld = vf.oldTime().oldTime();
            forAll(ddt0bf, patchI)
            {
                fvPatchField<Type>& pddt0 = ddt0bf[patchI];
                const fvPatchScalarField& prhoOld = rhoOld.boundaryField()[patchI];
                const fvPatchScalarField& prhoOldOld = rhoOldOld.boundaryField()[patchI];
                const fvPatchField<Type>& pvfOld = vfOld.boundaryField()[patchI];
                const fvPatchField<Type>& pvfOldOld = vfOldOld.boundaryField()[patchI];
                const labelUList& faceCells = pddt0.patch().faceCells();
                forAll(pddt0, i)
                {
                    label own = faceCells[i];
                    pddt0[i] = rDtCoef0[own] * (prhoOld[i] * pvfOld[i] - prhoOldOld[i] * pvfOldOld[i]);
                }
            }
        }

        forAll(source, cellI)
        {
            source[cellI] =
            (
                rDtCoef[cellI]*rho.oldTime().primitiveField()[cellI]
               *vf.oldTime().primitiveField()[cellI]
              + ocCoeff_[cellI]*ddt0.primitiveField()[cellI]
            )*mesh().V0()[cellI];
        }
    }
    else
    {
        if (evaluate(ddt0))
        {
            tmp<scalarField> trDtCoef0 = rDtCoef0_(ddt0);
            const scalarField& rDtCoef0 = trDtCoef0();

            forAll(ddt0.primitiveFieldRef(), cellI)
            {
                ddt0.primitiveFieldRef()[cellI] =
                    rDtCoef0[cellI]*
                    (
                        rho.oldTime().primitiveField()[cellI]
                       *vf.oldTime().primitiveField()[cellI]
                      - rho.oldTime().oldTime().primitiveField()[cellI]
                       *vf.oldTime().oldTime().primitiveField()[cellI]
                    )
                  - ocCoeff_[cellI]*ddt0.primitiveField()[cellI];
            }
        }

        forAll(source, cellI)
        {
            source[cellI] =
            (
                rDtCoef[cellI]*rho.oldTime().primitiveField()[cellI]
               *vf.oldTime().primitiveField()[cellI]
              + ocCoeff_[cellI]*ddt0.primitiveField()[cellI]
            )*V[cellI];
        }
    }

    return tfvm;
}


template<class Type>
tmp<fvMatrix<Type>>
momentumConsistentCrankNicolson<Type>::fvmDdt
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    DDt0Field<GeometricField<Type, fvPatchField, volMesh>>& ddt0 =
        ddt0_<GeometricField<Type, fvPatchField, volMesh>>
        (
            "ddt0(" + alpha.name() + ',' + rho.name() + ',' + vf.name() + ')',
            alpha.dimensions()*rho.dimensions()*vf.dimensions()
        );

    tmp<fvMatrix<Type>> tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            alpha.dimensions()*rho.dimensions()*vf.dimensions()*dimVol/dimTime
        )
    );
    fvMatrix<Type>& fvm = tfvm.ref();

    updateOcCoeff();
    tmp<scalarField> trDtCoef = rDtCoef_(ddt0);
    const scalarField& rDtCoef = trDtCoef();

    const scalarField& V = mesh().V();
    const scalarField& alphaField = alpha.primitiveField();
    const scalarField& rhoField = rho.primitiveField();
    scalarField& diag = fvm.diag();
    Field<Type>& source = fvm.source();

    forAll(diag, cellI)
    {
        diag[cellI] = rDtCoef[cellI]*alphaField[cellI]*rhoField[cellI]*V[cellI];
    }

    vf.oldTime().oldTime();
    alpha.oldTime().oldTime();
    rho.oldTime().oldTime();

    if (mesh().moving())
    {
        if (evaluate(ddt0))
        {
            tmp<scalarField> trDtCoef0 = rDtCoef0_(ddt0);
            const scalarField& rDtCoef0 = trDtCoef0();

            forAll(ddt0.primitiveFieldRef(), cellI)
            {
                ddt0.primitiveFieldRef()[cellI] =
                (
                    rDtCoef0[cellI]*
                    (
                        mesh().V0()[cellI]
                       *alpha.oldTime().primitiveField()[cellI]
                       *rho.oldTime().primitiveField()[cellI]
                       *vf.oldTime().primitiveField()[cellI]

                      - mesh().V00()[cellI]
                       *alpha.oldTime().oldTime().primitiveField()[cellI]
                       *rho.oldTime().oldTime().primitiveField()[cellI]
                       *vf.oldTime().oldTime().primitiveField()[cellI]
                    )
                  - mesh().V00()[cellI]
                   *ocCoeff_[cellI]*ddt0.primitiveField()[cellI]
                )/mesh().V0()[cellI];
            }

            typename GeometricField<Type, fvPatchField, volMesh>::Boundary& ddt0bf = ddt0.boundaryFieldRef();
            const volScalarField& alphaOld = alpha.oldTime();
            const volScalarField& alphaOldOld = alpha.oldTime().oldTime();
            const volScalarField& rhoOld = rho.oldTime();
            const volScalarField& rhoOldOld = rho.oldTime().oldTime();
            const GeometricField<Type, fvPatchField, volMesh>& vfOld = vf.oldTime();
            const GeometricField<Type, fvPatchField, volMesh>& vfOldOld = vf.oldTime().oldTime();
            forAll(ddt0bf, patchI)
            {
                fvPatchField<Type>& pddt0 = ddt0bf[patchI];
                const fvPatchScalarField& palphaOld = alphaOld.boundaryField()[patchI];
                const fvPatchScalarField& palphaOldOld = alphaOldOld.boundaryField()[patchI];
                const fvPatchScalarField& prhoOld = rhoOld.boundaryField()[patchI];
                const fvPatchScalarField& prhoOldOld = rhoOldOld.boundaryField()[patchI];
                const fvPatchField<Type>& pvfOld = vfOld.boundaryField()[patchI];
                const fvPatchField<Type>& pvfOldOld = vfOldOld.boundaryField()[patchI];
                const labelUList& faceCells = pddt0.patch().faceCells();
                forAll(pddt0, i)
                {
                    label own = faceCells[i];
                    pddt0[i] = rDtCoef0[own] *
                               (palphaOld[i] * prhoOld[i] * pvfOld[i] -
                                palphaOldOld[i] * prhoOldOld[i] * pvfOldOld[i]);
                }
            }
        }

        forAll(source, cellI)
        {
            source[cellI] =
            (
                rDtCoef[cellI]
               *alpha.oldTime().primitiveField()[cellI]
               *rho.oldTime().primitiveField()[cellI]
               *vf.oldTime().primitiveField()[cellI]
              + ocCoeff_[cellI]*ddt0.primitiveField()[cellI]
            )*mesh().V0()[cellI];
        }
    }
    else
    {
        if (evaluate(ddt0))
        {
            tmp<scalarField> trDtCoef0 = rDtCoef0_(ddt0);
            const scalarField& rDtCoef0 = trDtCoef0();

            forAll(ddt0.primitiveFieldRef(), cellI)
            {
                ddt0.primitiveFieldRef()[cellI] =
                    rDtCoef0[cellI]*
                    (
                        alpha.oldTime().primitiveField()[cellI]
                       *rho.oldTime().primitiveField()[cellI]
                       *vf.oldTime().primitiveField()[cellI]

                      - alpha.oldTime().oldTime().primitiveField()[cellI]
                       *rho.oldTime().oldTime().primitiveField()[cellI]
                       *vf.oldTime().oldTime().primitiveField()[cellI]
                    )
                  - ocCoeff_[cellI]*ddt0.primitiveField()[cellI];
            }
        }

        forAll(source, cellI)
        {
            source[cellI] =
            (
                rDtCoef[cellI]
               *alpha.oldTime().primitiveField()[cellI]
               *rho.oldTime().primitiveField()[cellI]
               *vf.oldTime().primitiveField()[cellI]
              + ocCoeff_[cellI]*ddt0.primitiveField()[cellI]
            )*V[cellI];
        }
    }

    return tfvm;
}


template<class Type>
tmp<typename momentumConsistentCrankNicolson<Type>::fluxFieldType>
momentumConsistentCrankNicolson<Type>::fvcDdtUfCorr
(
    const GeometricField<Type, fvPatchField, volMesh>& U,
    const GeometricField<Type, fvsPatchField, surfaceMesh>& Uf
)
{
    DDt0Field<GeometricField<Type, fvPatchField, volMesh>>& ddt0 =
        ddt0_<GeometricField<Type, fvPatchField, volMesh>>
        (
            "ddtCorrDdt0(" + U.name() + ')',
            U.dimensions()
        );

    DDt0Field<GeometricField<Type, fvsPatchField, surfaceMesh>>& dUfdt0 =
        ddt0_<GeometricField<Type, fvsPatchField, surfaceMesh>>
        (
            "ddtCorrDdt0(" + Uf.name() + ')',
            Uf.dimensions()
        );

    updateOcCoeff();
    tmp<scalarField> trDtCoef = rDtCoef_(ddt0);
    const scalarField& rDtCoef = trDtCoef();

    if (evaluate(ddt0))
    {
        tmp<scalarField> trDtCoef0 = rDtCoef0_(ddt0);
        const scalarField& rDtCoef0 = trDtCoef0();

        forAll(ddt0.primitiveFieldRef(), cellI)
        {
            ddt0.primitiveFieldRef()[cellI] =
                rDtCoef0[cellI]*
                (
                    U.oldTime().primitiveField()[cellI]
                  - U.oldTime().oldTime().primitiveField()[cellI]
                )
              - ocCoeff_[cellI]*ddt0.primitiveField()[cellI];
        }
    }

    if (evaluate(dUfdt0))
    {
        tmp<scalarField> trDtCoef0Uf = rDtCoef0_(dUfdt0);
        const scalarField& rDtCoef0Uf = trDtCoef0Uf();

        forAll(dUfdt0.primitiveFieldRef(), faceI)
        {
            label own = mesh().owner()[faceI];
            dUfdt0.primitiveFieldRef()[faceI] =
                rDtCoef0Uf[own]*
                (
                    Uf.oldTime().primitiveField()[faceI]
                  - Uf.oldTime().oldTime().primitiveField()[faceI]
                )
              - ocCoeff_[own]*dUfdt0.primitiveField()[faceI];
        }
    }

    tmp<GeometricField<Type, fvPatchField, volMesh>> trDtU
    (
        new GeometricField<Type, fvPatchField, volMesh>
        (
            IOobject
            (
                "rDtU",
                mesh().time().timeName(),
                mesh()
            ),
            mesh(),
            U.dimensions()/dimTime
        )
    );
    GeometricField<Type, fvPatchField, volMesh>& rDtU = trDtU.ref();

    forAll(rDtU.primitiveFieldRef(), cellI)
    {
        rDtU.primitiveFieldRef()[cellI] =
            rDtCoef[cellI]*U.oldTime().primitiveField()[cellI]
          + ocCoeff_[cellI]*ddt0.primitiveField()[cellI];
    }

    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> trDtUf
    (
        new GeometricField<Type, fvsPatchField, surfaceMesh>
        (
            IOobject
            (
                "rDtUf",
                mesh().time().timeName(),
                mesh()
            ),
            mesh(),
            Uf.dimensions()/dimTime
        )
    );
    GeometricField<Type, fvsPatchField, surfaceMesh>& rDtUf = trDtUf.ref();

    forAll(rDtUf.primitiveFieldRef(), faceI)
    {
        label own = mesh().owner()[faceI];
        rDtUf.primitiveFieldRef()[faceI] =
            rDtCoef[own]*Uf.oldTime().primitiveField()[faceI]
          + ocCoeff_[own]*dUfdt0.primitiveField()[faceI];
    }

    return tmp<fluxFieldType>
    (
        new fluxFieldType
        (
            IOobject
            (
                "ddtCorr(" + U.name() + ',' + Uf.name() + ')',
                mesh().time().timeName(),
                mesh().thisDb()
            ),
            this->fvcDdtPhiCoeff(U.oldTime(), mesh().Sf() & Uf.oldTime())
           *(
                mesh().Sf()
              & (
                    rDtUf
                  - fvc::interpolate(rDtU)
                )
            )
        )
    );
}


template<class Type>
tmp<typename momentumConsistentCrankNicolson<Type>::fluxFieldType>
momentumConsistentCrankNicolson<Type>::fvcDdtPhiCorr
(
    const GeometricField<Type, fvPatchField, volMesh>& U,
    const fluxFieldType& phi
)
{
    DDt0Field<GeometricField<Type, fvPatchField, volMesh>>& ddt0 =
        ddt0_<GeometricField<Type, fvPatchField, volMesh>>
        (
            "ddtCorrDdt0(" + U.name() + ')',
            U.dimensions()
        );

    DDt0Field<fluxFieldType>& dphidt0 =
        ddt0_<fluxFieldType>
        (
            "ddtCorrDdt0(" + phi.name() + ')',
            phi.dimensions()
        );

    updateOcCoeff();
    tmp<scalarField> trDtCoef = rDtCoef_(ddt0);
    const scalarField& rDtCoef = trDtCoef();

    if (evaluate(ddt0))
    {
        tmp<scalarField> trDtCoef0 = rDtCoef0_(ddt0);
        const scalarField& rDtCoef0 = trDtCoef0();

        forAll(ddt0.primitiveFieldRef(), cellI)
        {
            ddt0.primitiveFieldRef()[cellI] =
                rDtCoef0[cellI]*
                (
                    U.oldTime().primitiveField()[cellI]
                  - U.oldTime().oldTime().primitiveField()[cellI]
                )
              - ocCoeff_[cellI]*ddt0.primitiveField()[cellI];
        }
    }

    if (evaluate(dphidt0))
    {
        tmp<scalarField> trDtCoef0phi = rDtCoef0_(dphidt0);
        const scalarField& rDtCoef0phi = trDtCoef0phi();

        forAll(dphidt0.primitiveFieldRef(), faceI)
        {
            label own = mesh().owner()[faceI];
            dphidt0.primitiveFieldRef()[faceI] =
                rDtCoef0phi[own]*
                (
                    phi.oldTime().primitiveField()[faceI]
                  - phi.oldTime().oldTime().primitiveField()[faceI]
                )
              - ocCoeff_[own]*dphidt0.primitiveField()[faceI];
        }
    }

    tmp<GeometricField<Type, fvPatchField, volMesh>> trDtU
    (
        new GeometricField<Type, fvPatchField, volMesh>
        (
            IOobject
            (
                "rDtU",
                mesh().time().timeName(),
                mesh()
            ),
            mesh(),
            U.dimensions()/dimTime
        )
    );
    GeometricField<Type, fvPatchField, volMesh>& rDtU = trDtU.ref();

    forAll(rDtU.primitiveFieldRef(), cellI)
    {
        rDtU.primitiveFieldRef()[cellI] =
            rDtCoef[cellI]*U.oldTime().primitiveField()[cellI]
          + ocCoeff_[cellI]*ddt0.primitiveField()[cellI];
    }

    tmp<fluxFieldType> trDtPhi
    (
        new fluxFieldType
        (
            IOobject
            (
                "rDtPhi",
                mesh().time().timeName(),
                mesh()
            ),
            mesh(),
            phi.dimensions()/dimTime
        )
    );
    fluxFieldType& rDtPhi = trDtPhi.ref();

    forAll(rDtPhi.primitiveFieldRef(), faceI)
    {
        label own = mesh().owner()[faceI];
        rDtPhi.primitiveFieldRef()[faceI] =
            rDtCoef[own]*phi.oldTime().primitiveField()[faceI]
          + ocCoeff_[own]*dphidt0.primitiveField()[faceI];
    }

    return tmp<fluxFieldType>
    (
        new fluxFieldType
        (
            IOobject
            (
                "ddtCorr(" + U.name() + ',' + phi.name() + ')',
                mesh().time().timeName(),
                mesh().thisDb()
            ),
            this->fvcDdtPhiCoeff(U.oldTime(), phi.oldTime())
           *(
                rDtPhi
              - fvc::dotInterpolate(mesh().Sf(), rDtU)
            )
        )
    );
}


template<class Type>
tmp<typename momentumConsistentCrankNicolson<Type>::fluxFieldType>
momentumConsistentCrankNicolson<Type>::fvcDdtUfCorr
(
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& U,
    const GeometricField<Type, fvsPatchField, surfaceMesh>& Uf
)
{
    if
    (
        U.dimensions() == dimVelocity
     && Uf.dimensions() == rho.dimensions()*dimVelocity
    )
    {
        DDt0Field<GeometricField<Type, fvPatchField, volMesh>>& ddt0 =
            ddt0_<GeometricField<Type, fvPatchField, volMesh>>
            (
                "ddtCorrDdt0(" + rho.name() + ',' + U.name() + ')',
                rho.dimensions()*U.dimensions()
            );

        DDt0Field<GeometricField<Type, fvsPatchField, surfaceMesh>>& dUfdt0 =
            ddt0_<GeometricField<Type, fvsPatchField, surfaceMesh>>
            (
                "ddtCorrDdt0(" + Uf.name() + ')',
                Uf.dimensions()
            );

        updateOcCoeff();
        tmp<scalarField> trDtCoef = rDtCoef_(ddt0);
        const scalarField& rDtCoef = trDtCoef();

        GeometricField<Type, fvPatchField, volMesh> rhoU0
        (
            rho.oldTime()*U.oldTime()
        );

        if (evaluate(ddt0))
        {
            tmp<scalarField> trDtCoef0 = rDtCoef0_(ddt0);
            const scalarField& rDtCoef0 = trDtCoef0();

            forAll(ddt0.primitiveFieldRef(), cellI)
            {
                ddt0.primitiveFieldRef()[cellI] =
                    rDtCoef0[cellI]*
                    (
                        rhoU0.primitiveField()[cellI]
                      - rho.oldTime().oldTime().primitiveField()[cellI]
                       *U.oldTime().oldTime().primitiveField()[cellI]
                    )
                  - ocCoeff_[cellI]*ddt0.primitiveField()[cellI];
            }
        }

        if (evaluate(dUfdt0))
        {
            tmp<scalarField> trDtCoef0Uf = rDtCoef0_(dUfdt0);
            const scalarField& rDtCoef0Uf = trDtCoef0Uf();

            forAll(dUfdt0.primitiveFieldRef(), faceI)
            {
                label own = mesh().owner()[faceI];
                dUfdt0.primitiveFieldRef()[faceI] =
                    rDtCoef0Uf[own]*
                    (
                        Uf.oldTime().primitiveField()[faceI]
                      - Uf.oldTime().oldTime().primitiveField()[faceI]
                    )
                  - ocCoeff_[own]*dUfdt0.primitiveField()[faceI];
            }
        }

        tmp<GeometricField<Type, fvPatchField, volMesh>> trDtrhoU
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                IOobject
                (
                    "rDtrhoU",
                    mesh().time().timeName(),
                    mesh()
                ),
                mesh(),
                rho.dimensions()*U.dimensions()/dimTime
            )
        );
        GeometricField<Type, fvPatchField, volMesh>& rDtrhoU = trDtrhoU.ref();

        forAll(rDtrhoU.primitiveFieldRef(), cellI)
        {
            rDtrhoU.primitiveFieldRef()[cellI] =
                rDtCoef[cellI]*rhoU0.primitiveField()[cellI]
              + ocCoeff_[cellI]*ddt0.primitiveField()[cellI];
        }

        tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> trDtUf
        (
            new GeometricField<Type, fvsPatchField, surfaceMesh>
            (
                IOobject
                (
                    "rDtUf",
                    mesh().time().timeName(),
                    mesh()
                ),
                mesh(),
                Uf.dimensions()/dimTime
            )
        );
        GeometricField<Type, fvsPatchField, surfaceMesh>& rDtUf = trDtUf.ref();

        forAll(rDtUf.primitiveFieldRef(), faceI)
        {
            label own = mesh().owner()[faceI];
            rDtUf.primitiveFieldRef()[faceI] =
                rDtCoef[own]*Uf.oldTime().primitiveField()[faceI]
              + ocCoeff_[own]*dUfdt0.primitiveField()[faceI];
        }

        tmp<fluxFieldType> ddtCorr
        (
            new fluxFieldType
            (
                IOobject
                (
                    "ddtCorr("
                  + rho.name() + ',' + U.name() + ',' + Uf.name() + ')',
                    mesh().time().timeName(),
                    mesh().thisDb()
                ),
                this->fvcDdtPhiCoeff
                (
                    rhoU0,
                    mesh().Sf() & Uf.oldTime(),
                    rho.oldTime()
                )
               *(
                    mesh().Sf()
                  & (
                        rDtUf
                      - fvc::interpolate(rDtrhoU)
                    )
                )
            )
        );

        return ddtCorr;
    }
    else
    {
        FatalErrorInFunction
            << "dimensions of Uf are not correct"
            << abort(FatalError);

        return fluxFieldType::null();
    }
}


template<class Type>
tmp<typename momentumConsistentCrankNicolson<Type>::fluxFieldType>
momentumConsistentCrankNicolson<Type>::fvcDdtPhiCorr
(
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& U,
    const fluxFieldType& phi
)
{
    if
    (
        U.dimensions() == dimVelocity
     && phi.dimensions() == rho.dimensions()*dimVelocity*dimArea
    )
    {
        DDt0Field<GeometricField<Type, fvPatchField, volMesh>>& ddt0 =
            ddt0_<GeometricField<Type, fvPatchField, volMesh>>
            (
                "ddtCorrDdt0(" + rho.name() + ',' + U.name() + ')',
                rho.dimensions()*U.dimensions()
            );

        DDt0Field<fluxFieldType>& dphidt0 =
            ddt0_<fluxFieldType>
            (
                "ddtCorrDdt0(" + phi.name() + ')',
                phi.dimensions()
            );

        updateOcCoeff();
        tmp<scalarField> trDtCoef = rDtCoef_(ddt0);
        const scalarField& rDtCoef = trDtCoef();

        GeometricField<Type, fvPatchField, volMesh> rhoU0
        (
            rho.oldTime()*U.oldTime()
        );

        if (evaluate(ddt0))
        {
            tmp<scalarField> trDtCoef0 = rDtCoef0_(ddt0);
            const scalarField& rDtCoef0 = trDtCoef0();

            forAll(ddt0.primitiveFieldRef(), cellI)
            {
                ddt0.primitiveFieldRef()[cellI] =
                    rDtCoef0[cellI]*
                    (
                        rhoU0.primitiveField()[cellI]
                      - rho.oldTime().oldTime().primitiveField()[cellI]
                       *U.oldTime().oldTime().primitiveField()[cellI]
                    )
                  - ocCoeff_[cellI]*ddt0.primitiveField()[cellI];
            }
        }

        if (evaluate(dphidt0))
        {
            tmp<scalarField> trDtCoef0phi = rDtCoef0_(dphidt0);
            const scalarField& rDtCoef0phi = trDtCoef0phi();

            forAll(dphidt0.primitiveFieldRef(), faceI)
            {
                label own = mesh().owner()[faceI];
                dphidt0.primitiveFieldRef()[faceI] =
                    rDtCoef0phi[own]*
                    (
                        phi.oldTime().primitiveField()[faceI]
                      - phi.oldTime().oldTime().primitiveField()[faceI]
                    )
                  - ocCoeff_[own]*dphidt0.primitiveField()[faceI];
            }
        }

        tmp<GeometricField<Type, fvPatchField, volMesh>> trDtrhoU
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                IOobject
                (
                    "rDtrhoU",
                    mesh().time().timeName(),
                    mesh()
                ),
                mesh(),
                rho.dimensions()*U.dimensions()/dimTime
            )
        );
        GeometricField<Type, fvPatchField, volMesh>& rDtrhoU = trDtrhoU.ref();

        forAll(rDtrhoU.primitiveFieldRef(), cellI)
        {
            rDtrhoU.primitiveFieldRef()[cellI] =
                rDtCoef[cellI]*rhoU0.primitiveField()[cellI]
              + ocCoeff_[cellI]*ddt0.primitiveField()[cellI];
        }

        tmp<fluxFieldType> trDtPhi
        (
            new fluxFieldType
            (
                IOobject
                (
                    "rDtPhi",
                    mesh().time().timeName(),
                    mesh()
                ),
                mesh(),
                phi.dimensions()/dimTime
            )
        );
        fluxFieldType& rDtPhi = trDtPhi.ref();

        forAll(rDtPhi.primitiveFieldRef(), faceI)
        {
            label own = mesh().owner()[faceI];
            rDtPhi.primitiveFieldRef()[faceI] =
                rDtCoef[own]*phi.oldTime().primitiveField()[faceI]
              + ocCoeff_[own]*dphidt0.primitiveField()[faceI];
        }

        tmp<fluxFieldType> ddtCorr
        (
            new fluxFieldType
            (
                IOobject
                (
                    "ddtCorr("
                  + rho.name() + ',' + U.name() + ',' + phi.name() + ')',
                    mesh().time().timeName(),
                    mesh().thisDb()
                ),
                this->fvcDdtPhiCoeff(rhoU0, phi.oldTime(), rho.oldTime())
               *(
                    rDtPhi
                  - fvc::dotInterpolate(mesh().Sf(), rDtrhoU)
                )
            )
        );

        return ddtCorr;
    }
    else
    {
        FatalErrorInFunction
            << "dimensions of phi are not correct"
            << abort(FatalError);

        return fluxFieldType::null();
    }
}


template<class Type>
tmp<surfaceScalarField> momentumConsistentCrankNicolson<Type>::meshPhi
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    DDt0Field<surfaceScalarField>& meshPhi0 = ddt0_<surfaceScalarField>
    (
        "meshPhiCN_0",
        dimVolume
    );

    updateOcCoeff();
    tmp<scalarField> trCoef = coef_(meshPhi0);
    tmp<scalarField> trCoef0 = coef0_(meshPhi0);

    if (evaluate(meshPhi0))
    {
        const scalarField& coef0 = trCoef0();

        forAll(meshPhi0.primitiveFieldRef(), faceI)
        {
            label own = mesh().owner()[faceI];
            meshPhi0.primitiveFieldRef()[faceI] =
                coef0[own]*mesh().phi().oldTime().primitiveField()[faceI]
              - ocCoeff_[own]*meshPhi0.primitiveField()[faceI];
        }
    }

    tmp<surfaceScalarField> tmeshPhi
    (
        new surfaceScalarField
        (
            IOobject
            (
                mesh().phi().name(),
                mesh().time().timeName(),
                mesh().thisDb(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            mesh(),
            dimVolume/dimTime
        )
    );

    surfaceScalarField& meshPhi = tmeshPhi.ref();
    const scalarField& coef = trCoef();

    forAll(meshPhi.primitiveFieldRef(), faceI)
    {
        label own = mesh().owner()[faceI];
        meshPhi.primitiveFieldRef()[faceI] =
            coef[own]*mesh().phi().primitiveField()[faceI]
          - ocCoeff_[own]*meshPhi0.primitiveField()[faceI];
    }

    return tmeshPhi;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * Explicit Template Instantiation * * * * * * * * * //

#include "volFields.H"

namespace Foam
{
namespace fv
{
    // Instantiate the new fvcDdt overload for scalar fields
    template
    tmp<GeometricField<scalar, fvPatchField, volMesh>>
    momentumConsistentCrankNicolson<scalar>::fvcDdt
    (
        const volScalarField& alpha,
        const volScalarField& rho,
        const GeometricField<scalar, fvPatchField, volMesh>& vf
    );

    // Instantiate the new fvcDdt overload for vector fields
    template
    tmp<GeometricField<vector, fvPatchField, volMesh>>
    momentumConsistentCrankNicolson<vector>::fvcDdt
    (
        const volScalarField& alpha,
        const volScalarField& rho,
        const GeometricField<vector, fvPatchField, volMesh>& vf
    );
    
    // Instantiate the new fvmDdt overload for scalar fields
    template
    tmp<fvMatrix<scalar>>
    momentumConsistentCrankNicolson<scalar>::fvmDdt
    (
        const volScalarField& alpha,
        const volScalarField& rho,
        const GeometricField<scalar, fvPatchField, volMesh>& vf
    );

    // Instantiate the new fvmDdt overload for vector fields
    template
    tmp<fvMatrix<vector>>
    momentumConsistentCrankNicolson<vector>::fvmDdt
    (
        const volScalarField& alpha,
        const volScalarField& rho,
        const GeometricField<vector, fvPatchField, volMesh>& vf
    );
}
}

// ************************************************************************* //
