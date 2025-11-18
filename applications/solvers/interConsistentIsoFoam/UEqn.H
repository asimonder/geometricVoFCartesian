MRF.correctBoundaryVelocity(U);

tmp<volVectorField> tgradAlpha = fvc::grad(alpha1);
const volVectorField& gradAlpha = tgradAlpha();

forAll(intCoeff, cellI)
{
    if (mag(gradAlpha[cellI]) > interfaceThreshold)
    {
        intCoeff[cellI] = 0.0; 
    }
    else
    {
        intCoeff[cellI] = 1.0; 
    }
}

//surfScalarField impCoeffs = fvc::interpolate(1.0 - (1.0 - beta)*intCoeff);
//surfScalarField expCoeffs = fvc::interpolate(1.0 - (1.0 - beta)*intCoeff);
tmp<surfaceScalarField> timpCoeffs = fvc::interpolate(1.0 - (1.0 - beta)*intCoeff);
const surfaceScalarField& impCoeffs = timpCoeffs();
tmp<surfaceScalarField> texpCoeffs = fvc::interpolate((1.0 - beta)*intCoeff);
const surfaceScalarField& expCoeffs = texpCoeffs();


fvVectorMatrix UEqn
(
    fvm::ddt(rho, U) + MRF.DDt(rho, U)
      +  fvm::div(impCoeffs*rhoPhi, U)
      - fvm::laplacian(impCoeffs*mixture.muf(), U)
  ==
    fvOptions(rho, U)
    - fvc::div(expCoeffs*rhoPhi.oldTime(), U.oldTime())
    + fvc::laplacian(expCoeffs*mufOld, U.oldTime())
);


UEqn.relax();

fvOptions.constrain(UEqn);

if (pimple.momentumPredictor())
  {
    solve
      (
       UEqn
       ==
       fvc::reconstruct
       (
	(
	 interfaceForce.surfaceTensionForce()
	 - ghf*fvc::snGrad(rho)
	 - fvc::snGrad(p_rgh)
	 ) * mesh.magSf()
	)
       );
    
    fvOptions.correct(U);
  }

