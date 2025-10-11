/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Computes the reconstructed distance function (surfaceDistance) from a volume-fraction
    field using geometricVoF reconstruction.

\*---------------------------------------------------------------------------*/

#include "surfaceDistance.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(surfaceDistance, 0);
    addToRunTimeSelectionTable(functionObject, surfaceDistance, dictionary);
}
}


// * * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * //
bool Foam::functionObjects::surfaceDistance::calc()
{
    if (!obr().foundObject<volScalarField>(alphaName_))
    {
        WarningInFunction
            << "alpha field not found: " << alphaName_ << nl;
        return false;
    }

    const volScalarField& alpha =
        obr().lookupObject<volScalarField>(alphaName_);

    const fvMesh& mesh = alpha.mesh();
    const ijkZone& IJK = ijkZone::New(mesh);

    const label Nx = IJK.Nx();
    const label Ny = IJK.Ny();
    const scalar dz = IJK.dz();

    auto ij = [&](const label i, const label j) -> label { return j*Nx + i; };

    // integrate alpha vertically for each column
    Field<scalar> eta(Nx*Ny, 0.0);

    forAll(alpha, celli)
    {
        const Vector<label> g = IJK.ijk3(celli);
        const label i = g.x(), j = g.y();

        if (i < 0 || i >= Nx || j < 0 || j >= Ny)
            continue;

        eta[ij(i,j)] += alpha[celli] * dz;
    }

    // sum over ranks
    Pstream::reduce(eta, sumOp<scalar>());

    // construct the resulting volScalarField
    tmp<volScalarField> tPhi
    (
        new volScalarField
        (
            IOobject
            (
                resultName_,                // from dict, e.g. "surfaceDistance"
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("zero", dimLength, 0.0)
        )
    );

    volScalarField& phi = tPhi.ref();

    // assign distances: η(i,j) - z_c
    forAll(phi, celli)
    {
        const Vector<label> g = IJK.ijk3(celli);
        const label i = g.x(), j = g.y();
        if (i < 0 || i >= Nx || j < 0 || j >= Ny)
            continue;

        phi[celli] = eta[ij(i,j)] - mesh.C()[celli].z();
    }

    phi.correctBoundaryConditions();

    // store and return true on success
    return store(resultName_, tPhi);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::surfaceDistance::surfaceDistance
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
  //fvMeshFunctionObject(name, runTime, dict)
    fieldExpression(name, runTime, dict, "alpha.water")
  //  fieldExpression(name, runTime, dict)
    // alphaName_(dict.lookupOrDefault<word>("alpha", "alpha.water"))
{
  setResultName(typeName,"");
  //resultName_ = typeName; 
}


// ************************************************************************* //
