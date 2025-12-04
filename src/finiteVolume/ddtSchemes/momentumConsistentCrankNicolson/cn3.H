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

//dimensionedScalar rDtCoef = rDtCoef_(ddt0);
updateOcCoeff();
//tmp<scalarField> trDtCoef = rDtCoef_(ddt0)().primitiveField();
//  const scalarField& rDtCoef = trDtCoef();

scalar phiOldAv = gSum(mag(phi.oldTime().primitiveField())*mesh().magSf())/gSum(mesh().magSf());
Info<< "phi.oldTime() average: " << phiOldAv << endl;
scalar phiOldOldAv = gSum(mag(phi.oldTime().oldTime().primitiveField())*mesh().magSf())/gSum(mesh().magSf());
Info<< "phi.oldTime().oldTime() average: " << phiOldOldAv << endl;

if (evaluate(ddt0))
{
    ddt0 =
      rDtCoef0_(ddt0)*(U.oldTime() - U.oldTime().oldTime())
      - offCentre_(ddt0());

    /*    scalar offCenterAv = gSum(mag(offCentre_(ddt0())().primitiveField())*mesh().V())/gSum(mesh().V());
    Info<< "offCenter average: " << offCenterAv << endl;
    scalar rdtAv = gSum(mag(rDtCoef0_(ddt0)().primitiveField())*mesh().V())/gSum(mesh().V());
    Info<< "rDtCoef0_(ddt0) average: " << rdtAv << endl;
    scalar ddt0MagAv = gSum(mag(ddt0().primitiveField())*mesh().V())/gSum(mesh().V());
    Info<< "Mag ddt0 average: " << ddt0MagAv << endl;*/
}

if (evaluate(dphidt0))
{
    // dphidt0 =
    //  rDtCoef0_(dphidt0)*(phi.oldTime() - phi.oldTime().oldTime())
    //  - offCentre_(dphidt0());
    dphidt0 =
      fvc::interpolate(rDtCoef0_(ddt0))*(phi.oldTime() - phi.oldTime().oldTime())
      //- offCentre_(dphidt0());
      -fvc::interpolate(ocCoeff_)*dphidt0();

    /*    scalar interpRdtAv = gSum(mag(fvc::interpolate(rDtCoef0_(ddt0))().primitiveField())*mesh().magSf())/gSum(mesh().magSf());
    Info<< "interpolated rDtCoef0_(ddt0) average: " << interpRdtAv << endl;
    scalar interpOcAv = gSum(mag(fvc::interpolate(ocCoeff_)().primitiveField())*mesh().magSf())/gSum(mesh().magSf());
    Info<< "interpolated ocCoeff_ average: " << interpOcAv << endl;
    scalar phiDiffAv = gSum(mag((phi.oldTime() - phi.oldTime().oldTime())().primitiveField())*mesh().magSf())/gSum(mesh().magSf());
    Info<< "phi difference average: " << phiDiffAv << endl;
    scalar dphidt0Mag = gSum(mag(dphidt0().primitiveField())*mesh().magSf())/gSum(mesh().magSf());
    Info<< "dphidt0 mag (global): " << dphidt0Mag << endl;*/
    }

dimensionedScalar rDtCoefScalar("rDt", dimensionSet(0, 0, -1, 0, 0, 0, 0), 1.0);
surfaceScalarField rDtCoeff(fvc::interpolate(rDtCoef_(ddt0)));
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
        0.0*this->fvcDdtPhiCoeff(U.oldTime(), phi.oldTime())
       *(
         //(rDtCoef*phi.oldTime() + offCentre_(dphidt0()))
	  (fvc::interpolate(rDtCoef_(ddt0))*phi.oldTime() + fvc::interpolate(ocCoeff_)*dphidt0())
	 - fvc::dotInterpolate
            (
                mesh().Sf(),
                rDtCoeff*U.oldTime() + offCentre_(ddt0())
	   )
        )
    )
);

//scalar ddtCorrMag = gSum(mag(ddtCorr().primitiveField())*mesh().magSf())/gSum(mesh().magSf());
//Info<< "ddtCorr mag (global): " << ddtCorrMag << endl;

return ddtCorr;
