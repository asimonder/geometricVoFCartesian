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

tmp<volScalarField> tOffCenterCoeff = beta * intCoeff;
const volScalarField& offCenterCoeff = tOffCenterCoeff();

if (runTime.timeIndex() > runTime.startTimeIndex() + 1)
  {
    ddt0 = ((1.0 + offCenterCoeff) / dt0) * deltaRhoU - (offCenterCoeff * ddt0_old_stored);

    //ddt0 = rDtCoef0_(ddt0)*
    //  (
    //  rho.oldTime()*vf.oldTime()
    //  - rho.oldTime().oldTime()*vf.oldTime().oldTime()
    //  ) - offCentre_(ddt0());
    
  }
 else
   {
     ddt0 = deltaRhoU / dt0;
   }

fvVectorMatrix UEqn
(
    (1.0 + offCenterCoeff) * fvm::ddt(rho, U)
  + (1.0 + offCenterCoeff) * MRF.DDt(rho, U)
  + fvm::div(rhoPhi, U)
  - fvm::laplacian(mixture.muf(), U)
 ==
    fvOptions(rho, U)
  + offCenterCoeff * ddt0
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

