// self-contained interface Courant number (cells with mixed alpha)
scalar maxAlphaCo
(
    runTime.controlDict().get<scalar>("maxAlphaCo")
);

scalar alphaCoNum = 0.0;

if (mesh.nInternalFaces())
{
    scalarField sumPhi
    (
        fvc::surfaceSum(mag(phi))().primitiveField()
    );
    forAll(sumPhi, celli)
    {
        const scalar a = alpha1[celli];
        if (a < 1e-6 || a > 1 - 1e-6)
        {
            sumPhi[celli] = 0;
        }
    }
    alphaCoNum =
        0.5*gMax(sumPhi/mesh.V().field())*runTime.deltaTValue();
}

Info<< "Interface Courant Number max: " << alphaCoNum << endl;
