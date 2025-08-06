/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2021-2025 Asim Onder
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "ijkZone.H"
#include "exprOps.H"
#include "processorPolyPatch.H"
#include "processorCyclicPolyPatch.H"
#include "PstreamReduceOps.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(ijkZone, 0);
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
// This is the implementation of your new static function.
Foam::ijkZone& Foam::ijkZone::New(const fvMesh& mesh)
{
    // Check the mesh's object registry to see if an ijkZone already exists.
    if (!mesh.thisDb().foundObject<ijkZone>(ijkZone::typeName))
    {
        // If not, create one...
        Info<< "Creating new ijkZone object and adding to registry." << endl;
        ijkZone* ijkZonePtr = new ijkZone(mesh);
        // ...and store it in the registry.
        ijkZonePtr->store();
    }

    // Return a reference to the object now guaranteed to be in the registry.
    return const_cast<ijkZone&>(mesh.lookupObject<ijkZone>(ijkZone::typeName));
}


Foam::ijkZone::ijkZone(const fvMesh& mesh):
  regIOobject(IOobject
	      (
	       typeName,
	       mesh.time().system(), 
	       mesh,
	       IOobject::NO_READ,
	       IOobject::NO_WRITE,
	       false
	       )),
  mesh_(mesh),
  globalNumbering_(mesh_.nCells()+mesh_.nBoundaryFaces()),
  C_(mesh.C()),
  isEmpty_(false,false,false),
  isCyclic_(false,false,false)
{

  Info<<"Building the ijkZone..."<<endl;
  IOdictionary ijkDict
    (
     IOobject
     (
      "ijkDict",
      mesh_.time().system(),
      mesh_.time().db(),
      IOobject::MUST_READ,
      IOobject::NO_WRITE
      )
     );
  
  word cellSetName=ijkDict.get<word>("cellSet");
  bool isDomain=true;
  if (cellSetName!="domain")
    isDomain=false;

  symXIn_=ijkDict.lookupOrDefault<bool>("symXIn",false);
  symXOut_=ijkDict.lookupOrDefault<bool>("symXOut",false);
  symYIn_=ijkDict.lookupOrDefault<bool>("symYIn",false);
  symYOut_=ijkDict.lookupOrDefault<bool>("symYOut",false);
  symZIn_=ijkDict.lookupOrDefault<bool>("symZIn",false);
  symZOut_=ijkDict.lookupOrDefault<bool>("symZOut",false);
  isCyclic_.x()=ijkDict.lookupOrDefault<bool>("cyclicX",false);
  isCyclic_.y()=ijkDict.lookupOrDefault<bool>("cyclicY",false);
  isCyclic_.z()=ijkDict.lookupOrDefault<bool>("cyclicZ",false);

  Info<<"symXIn="<<symXIn_<<endl;
  
  const  surfaceVectorField&  Cf = mesh_.Cf();

  const cellList& cells=mesh_.cells();
  labelList cellIndices;
  if (!isDomain)
    {
      cellSet cellsToWrite(mesh_, cellSetName);    
      cellIndices=cellsToWrite.toc();
      Info<<"Found the cellSet "<<cellSetName<<" for construction of ijkZone..."<<endl;
    }
  else
    {
      Info<<"Using the whole domain for the construction of ijkZone..."<<endl;
      cellIndices.setSize(C_.size(),0.0);
      forAll(C_,celli)
	cellIndices[celli]=celli;
    }

  Info<<"Size of the interface domain of the master proc="<<cellIndices.size()<<endl;
  Pmin_=point::uniform(VGREAT);
  Pmax_=point::uniform(-VGREAT);
  dx_=0.0;dy_=0.0;dz_=0.0;

  if (cellIndices.size()>0)
    {
      const faceList& faces=mesh_.faces();
      const label iCell=cellIndices[0];
      for (int iF=0;iF<6;iF++)
       {	    
         label f1=cells[iCell][iF];
         label f2=cells[iCell].opposingFaceLabel(f1,faces);
         scalar dx=Cf[f1].x()-Cf[f2].x();
	 scalar dy=Cf[f1].y()-Cf[f2].y();
         scalar dz=Cf[f1].z()-Cf[f2].z();
         if (dx>dx_)
	   dx_=dx;
	 if (dy>dy_)
	   dy_=dy;
	 if (dz>dz_)
	   dz_=dz;		    
	}

      forAll(cellIndices,celli)
	{
	  const label iCell=cellIndices[celli];
	  const point cc = C_[iCell];
	  Pmin_=Foam::min(Pmin_, cc);
	  Pmax_=Foam::max(Pmax_, cc);
	}
    }
  
  if (Pstream::parRun())
    {
      Foam::reduce(Pmin_, minOp<point>());
      Foam::reduce(Pmax_, maxOp<point>());
      Foam::reduce(dx_, maxOp<scalar>());
      Foam::reduce(dy_, maxOp<scalar>());
      Foam::reduce(dz_, maxOp<scalar>());    
    }

  Info<< "dx= "<<dx_<< ",  dy= "<<dy_<< ",   dz= "<<dz_<<endl;

  Lx_=Pmax_.x()-Pmin_.x();
  Ly_=Pmax_.y()-Pmin_.y();
  Lz_=Pmax_.z()-Pmin_.z();
  Nx_=round(Lx_/dx_+1.0);
  Ny_=round(Ly_/dy_+1.0);
  Nz_=round(Lz_/dz_+1.0);
  

  Info<<"Pmin="<<Pmin_<<", Pmax="<<Pmax_<<endl;	   
  Info<<"Nx= "<<Nx_<<", Ny= "<<Ny_<<", Nz= "<<Nz_<<endl;

  if (Nx_==1)
    isEmpty_.x()=true;
  if (Ny_==1)
    isEmpty_.y()=true;
  if (Nz_==1)
    isEmpty_.z()=true;

  Info<<"is x periodic: " <<isCyclic_.x()<<endl;
  Info<<"is y periodic: " <<isCyclic_.y()<<endl;
  Info<<"is z periodic: " <<isCyclic_.z()<<endl;

  
 ///////////////////// REFACTORING THE PARALLEL COMMUNICATION ///////////////////////
#include "parallelCommunication.H"

}

//Foam::List<bool> Foam::ijkZone::markBoundaryCells(label nOffsetCells)
void Foam::ijkZone::markBoundaryCells(boolList& atBoundary,label nOffsetCells)
{
  //boolList atBoundary(mesh_.nCells(),false);
  forAll(C_,celli)
    {
      Vector<label> ijk=ijk3(celli);
      if ((!isCyclic_.x() and !isEmpty_.x() and ((ijk.x()+nOffsetCells>=Nx_ and !symXOut_) or (ijk.x()-nOffsetCells<=0 and !symXIn_))) or
	  (!isCyclic_.y() and !isEmpty_.y() and ((ijk.y()+nOffsetCells>=Ny_ and !symYOut_) or (ijk.y()-nOffsetCells<=0 and !symYIn_))) or
	  (!isCyclic_.z() and !isEmpty_.z() and ((ijk.z()+nOffsetCells>=Nz_ and !symZOut_) or (ijk.z()-nOffsetCells<=0 and !symZIn_))))
	atBoundary[celli]=true;
    }

  //return atBoundary;      

}


Foam::Vector<Foam::label> Foam::ijkZone::ijk3(label celli)
{  
  label i=round((C_[celli].x()-Pmin_.x())/dx_);
  label j=round((C_[celli].y()-Pmin_.y())/dy_);
  label k=round((C_[celli].z()-Pmin_.z())/dz_);

  if (isEmpty_.x())
    i=0;
  else if (isEmpty_.y())
    j=0;
  else if (isEmpty_.z())
    k=0;

  return Vector<label>(i,j,k);
}



// * * * * * * * * * * * * * * * * Selectors  * * * * * * * * * * * * * * //
/*Foam::ijkZone& Foam::ijkZone::New(const fvMesh& mesh)
{
    bool found = mesh.thisDb().foundObject<ijkZone>
    (
        ijkZone::typeName
    ); 
    ijkZone* ptr = nullptr;

    if(found)
    {
        ptr = mesh.thisDb().getObjectPtr<ijkZone>
        (
            ijkZone::typeName
        );

        return *ptr;
    }

    ijkZone* objectPtr = new ijkZone(mesh);

    regIOobject::store(static_cast<ijkZone*>(objectPtr));

    return *objectPtr;
    }*/
 
bool Foam::ijkZone::writeData(Ostream& os) const
{
    // This object has no data to write, but the function must exist.
    // We could write the globalIds_ map here if we wanted to, but it's not necessary.
    return os.good();
}

// ************************************************************************* //
