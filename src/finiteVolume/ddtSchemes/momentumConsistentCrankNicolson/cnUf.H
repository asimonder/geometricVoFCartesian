// Dynamic-mesh variant of the transient flux correction (ddtCorr) with
// the spatially varying off-centering coefficient. Mirrors cn2.H with
// Uf in place of phi. AMR and moving-mesh solvers reach this path,
// because fvc::ddtCorr(U, phi, Uf) dispatches to the Uf form whenever
// the mesh is dynamic.

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

if (evaluate(ddt0))
{
    ddt0 =
        rDtCoef0_(ddt0)*(U.oldTime() - U.oldTime().oldTime())
      - offCentre_(ddt0());
}

if (evaluate(dUfdt0))
{
    dUfdt0 =
        fvc::interpolate(rDtCoef0_(ddt0))
       *(Uf.oldTime() - Uf.oldTime().oldTime())
      - fvc::interpolate(ocCoeff_)*dUfdt0();
}

tmp<fluxFieldType> ddtCorr
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
            (
                mesh().Sf()
              & (
                    fvc::interpolate(rDtCoef_(ddt0))*Uf.oldTime()
                  + fvc::interpolate(ocCoeff_)*dUfdt0()
                )
            )
          - fvc::dotInterpolate
            (
                mesh().Sf(),
                rDtCoef_(ddt0)*U.oldTime() + offCentre_(ddt0())
            )
        )
    )
);

return ddtCorr;
