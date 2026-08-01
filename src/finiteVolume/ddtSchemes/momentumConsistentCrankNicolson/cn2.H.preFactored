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
dphidt0.setOriented();

updateOcCoeff();

if (evaluate(ddt0))
{
    ddt0 =
      rDtCoef0_(ddt0)*(U.oldTime() - U.oldTime().oldTime())
      - offCentre_(ddt0());

}

if (evaluate(dphidt0))
{
    dphidt0 =
      fvc::interpolate(rDtCoef0_(ddt0))*(phi.oldTime() - phi.oldTime().oldTime())
      -fvc::interpolate(ocCoeff_)*dphidt0();

    }

dimensionedScalar rDtCoefScalar("rDt", dimensionSet(0, 0, -1, 0, 0, 0, 0), 1.0);
tmp<fluxFieldType> ddtCorr
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
	  (fvc::interpolate(rDtCoef_(ddt0))*phi.oldTime() + fvc::interpolate(ocCoeff_)*dphidt0())
	 - fvc::dotInterpolate
            (
                mesh().Sf(),
                rDtCoef_(ddt0)*U.oldTime() + offCentre_(ddt0())
	   )
        )
    )
);

return ddtCorr;
