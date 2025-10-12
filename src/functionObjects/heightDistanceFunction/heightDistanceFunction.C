/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Computes the reconstructed distance function (heightDistanceFunction) from a volume-fraction
    field using geometricVoF reconstruction.

\*---------------------------------------------------------------------------*/

#include "heightDistanceFunction.H"
#include "addToRunTimeSelectionTable.H"
#include "PstreamReduceOps.H"

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(heightDistanceFunction, 0);
    addToRunTimeSelectionTable(functionObject, heightDistanceFunction, dictionary);
}
}


// * * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * //
bool Foam::functionObjects::heightDistanceFunction::calc()
{
    if (!obr().foundObject<volScalarField>(alphaName_))
    {
        WarningInFunction
            << "alpha field not found: " << alphaName_ << nl;
        return false;
    }

    const volScalarField& alpha =
        obr().lookupObject<volScalarField>(alphaName_);

    auto ij = [&](const label i, const label j) -> label { return j*Nx_ + i; };

    // integrate alpha vertically for each column
    scalarField eta(Nx_*Ny_, 0.0);

    forAll(alpha, celli)
    {
        const Vector<label> g = IJK_.ijk3(celli);
        const label i = g.x(), j = g.y();

	if (mesh_.C()[celli].z()<IJK_.Pmin().z() || mesh_.C()[celli].z()>IJK_.Pmax().z())
	  continue;

        eta[ij(i,j)] += alpha[celli] * dz_;
    }

    forAll(eta, i)
      {
	reduce(eta[i], sumOp<scalar>());
      }

    eta_=eta;
    
    // construct the resulting volScalarField
    tmp<volScalarField> tPhi
    (
        new volScalarField
        (
            IOobject
            (
                resultName_,                // from dict, e.g. "heightDistanceFunction"
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("zero", dimLength, 0.0)
        )
    );

    volScalarField& phi = tPhi.ref();

    forAll(phi, celli)
    {
        const Vector<label> g = IJK_.ijk3(celli);
        const label i = g.x(), j = g.y();
        if (i < 0 || i >= Nx_ || j < 0 || j >= Ny_)
            continue;

        phi[celli] = mesh_.C()[celli].z()-(eta[ij(i,j)]+zMin_);
    }

    phi.correctBoundaryConditions();
    
    // store and return true on success
    return store(resultName_, tPhi);
}


bool Foam::functionObjects::heightDistanceFunction::write()
{
    // First, call parent write (writes the field if needed)
    bool ok = fieldExpression::write();
    
    // Then write interface data
    if (Pstream::master())
    {
        fileName outputDir = 
            obr().time().globalPath()/
            "postProcessing"/
            name()/
            obr().time().timeName();
        
        mkDir(outputDir);
        
        OFstream os(outputDir/"interface.dat");

	os << "# Nx " << IJK_.Nx() << " Ny " << IJK_.Ny() << nl;
        os << "# x y eta" << nl;
        
        auto ij = [&](const label i, const label j) -> label { return j*Nx_ + i; };
        
        for (label j = 0; j < Ny_; j++)
        {
            for (label i = 0; i < Nx_; i++)
            {
                os << i*dx_+dx/2. << " " << j*dy_+dy_/2. << " " << eta_[ij(i,j)] << nl;
            }
        }
    }
    
    return ok;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::heightDistanceFunction::heightDistanceFunction
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
  //fvMeshFunctionObject(name, runTime, dict)
  fieldExpression(name, runTime, dict, "alpha.water"),
  IJK_(Foam::ijkZone::New(mesh_)),
  alphaName_(dict.lookupOrDefault<word>("alpha", "alpha.water"))
  //  meanWaterLevel_(dict.lookupOrDefault<scalar>("meanWaterLevel", 0.0))
  
{
  setResultName(typeName,"");
  //resultName_ = typeName; 

  Nx_ = IJK_.Nx();
  Ny_ = IJK_.Ny();
  dx_ = IJK_.dx();
  dy_ = IJK_.dy();
  dz_ = IJK_.dz();

  zMin_ = IJK_.Pmin().z()-(IJK_.dz()/2.);
    
}

// ************************************************************************* //
