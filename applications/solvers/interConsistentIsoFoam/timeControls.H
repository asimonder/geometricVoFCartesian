
// Look up the fvSolution dictionary from the database
//const IOdictionary& fvSolutionDict =
// mesh.db().lookupObject<IOdictionary>("fvSolution");
const IOdictionary& fvSolutionDict = mesh.lookupObject<IOdictionary>("fvSolution");

// Look up your custom sub-dictionary
const dictionary& numericalControlsDict =
  fvSolutionDict.subDict("timeControls");

// Read the 'beta' entry from the sub-dictionary
const scalar beta = readScalar(numericalControlsDict.lookup("beta"));
const scalar interfaceThreshold = readScalar(numericalControlsDict.lookup("interfaceThreshold"));

Info<< "Using custom temporal blending with beta = " << beta << endl;

volScalarField intCoeff
(
    IOobject
    (
        "intCoeff",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("intCoeff", dimless, 1.0) 
);
