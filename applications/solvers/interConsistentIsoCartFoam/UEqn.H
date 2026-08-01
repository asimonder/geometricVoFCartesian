MRF.correctBoundaryVelocity(U);

#include "configureInterfaceLayer.H"
#include "upwindDiagnostic.H"

// 3. Construct UEqn with BOTH terms
fvVectorMatrix UEqn
(
    fvm::ddt(rho, U) 
    + fvm::div(rhoPhiInter, U)  // Uses "upwind" from fvSchemes
    + fvm::div(rhoPhiBulk, U)    // Uses "limitedLinear" from fvSchemes
    + MRF.DDt(rho, U)
    + turbulence->divDevRhoReff(rho, U)
 ==
    fvOptions(rho, U)
);

/*fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U) + fvm::div(rhoPhi, U)
        + MRF.DDt(rho, U)
	+ turbulence->divDevRhoReff(rho, U)
     ==
        fvOptions(rho, U)
	);*/

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

