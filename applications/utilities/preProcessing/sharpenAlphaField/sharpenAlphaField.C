#include "fvCFD.H"

using namespace Foam;

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    Info<< "Reading sharpenAlphaDict\n" << endl;

    IOdictionary dict
    (
        IOobject
        (
            "sharpenAlphaDict",
	    runTime.globalPath()/"system",
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    const word fieldName = dict.lookupOrDefault<word>("alphaFieldName", "alpha.water");
    const scalar minAlpha = dict.lookupOrDefault<scalar>("alphaMinThreshold", -GREAT);
    const scalar maxAlpha = dict.lookupOrDefault<scalar>("alphaMaxThreshold", GREAT);

    Info<< "alphaFieldName           : " << fieldName << nl
        << "alphaMinThreshold        : " << minAlpha << nl
        << "alphaMaxThreshold        : " << maxAlpha << nl << endl;

    volScalarField alpha
    (
        IOobject
        (
            fieldName,
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    forAll(alpha, i)
    {
      if (alpha[i]<minAlpha)
	alpha[i]=0.0;
      else if (alpha[i]>maxAlpha)
	alpha[i]=1.0;
    }

    alpha.correctBoundaryConditions();
    alpha.write();

    Info<< "Finished sharpening alpha field " << fieldName << nl
        << "Cells below alphaMinThreshold set to 0.0, above alphaMaxThreshold set to 1.0." << endl;

    return 0;
}
