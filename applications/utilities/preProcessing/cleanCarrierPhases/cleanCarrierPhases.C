#include "fvCFD.H"

using namespace Foam;

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    Info<< "Reading cleanCarrierPhasesDict\n" << endl;

    IOdictionary dict
    (
        IOobject
        (
            "cleanCarrierPhasesDict",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    const word fieldName = dict.lookupOrDefault<word>("field", "alpha.water");
    const scalar zMinAlpha = dict.lookupOrDefault<scalar>("zMinAlpha", -GREAT);
    const scalar zMaxAlpha = dict.lookupOrDefault<scalar>("zMaxAlpha", GREAT);

    Info<< "Field           : " << fieldName << nl
        << "zMinAlpha       : " << zMinAlpha << nl
        << "zMaxAlpha       : " << zMaxAlpha << nl << endl;

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

    const scalarField& zc = mesh.C().component(vector::Z);

    forAll(alpha, i)
    {
        const scalar z = zc[i];
        if (z < zMinAlpha)
        {
            alpha[i] = 1.0;
        }
        else if (z > zMaxAlpha)
        {
            alpha[i] = 0.0;
        }
        // leave cells between zMinAlpha and zMaxAlpha unchanged
    }

    alpha.correctBoundaryConditions();
    alpha.write();

    Info<< "Finished cleaning " << fieldName << nl
        << "Cells below zMinAlpha set to 1.0, above zMaxAlpha set to 0.0." << endl;

    return 0;
}
