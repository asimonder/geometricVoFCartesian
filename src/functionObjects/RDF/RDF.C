/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Computes the reconstructed distance function (RDF) from a volume-fraction
    field using geometricVoF reconstruction.

\*---------------------------------------------------------------------------*/

#include "RDF.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(RDF, 0);
    addToRunTimeSelectionTable(functionObject, RDF, dictionary);
}
}


// * * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * //
bool Foam::functionObjects::RDF::calc()
{
  /*if (!obr().foundObject<volScalarField>(alphaName_))
    {
        WarningInFunction << "alpha field not found: " << alphaName_ << nl;
        return false;
	}*/

  //const volScalarField& alpha = obr().lookupObject<volScalarField>(alphaName_);
  

    reconstructionSchemes& recon =
      obr().lookupObjectRef<reconstructionSchemes>("reconstructionScheme");

  //recon.reconstruct(false);
  //recon.reconstruct(true);
  
    const volVectorField& interfaceCentres = recon.centre();
    const volVectorField& interfaceNormals = recon.normal();
    zoneDistribute& exchangeFields = zoneDistribute::New(mesh_);


    reconstructedDistanceFunction rdf(mesh_);
    boolList interfaceCells(mesh_.nCells(), false);
    forAll(interfaceCentres,cellI)
    {
        if (mag(interfaceNormals[cellI]) != 0)
        {
            interfaceCells[cellI] = true;
        }
        else
        {
            interfaceCells[cellI] = false;
        }
    }

    // rdf.markCellsNearSurf(recon.interfaceCell(),20);
    rdf.markCellsNearSurf(interfaceCells,20);
    const boolList& nextToInterface= rdf.nextToInterface();
    //boolList nextToInterface(mesh_.nCells(), true);
    //exchangeFields.updateStencil(nextToInterface);
    const Foam::volScalarField& phi=rdf.constructRDF(nextToInterface, interfaceCentres, interfaceNormals, exchangeFields, true);

    //    rdf.correctBoundaryConditions(); // update proc patches

    Info<< "Storing RDF field: min=" << gMin(phi) 
    << " max=" << gMax(phi) << endl;


    //    forAll(phi,iCell)
    // {
    //	if (phi[iCell]==0)
    //	  phi[iCell]=1000.0;
    // }
    
    //    return store(resultName_, phi);
    return store(resultName_, tmp<volScalarField>(new volScalarField(phi)));
}

//bool Foam::functionObjects::RDF::execute()
//{
    // Typically do nothing here, only act at write time
//  return true;
//}

//bool Foam::functionObjects::RDF::write()
//{
    // Call calc() to actually compute and write RDF at output times
//  return calc();
//}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::RDF::RDF
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
