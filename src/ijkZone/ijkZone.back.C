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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(ijkZone, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ijkZone::ijkZone(const fvMesh& mesh):
  //MeshObject<fvMesh, Foam::TopologicalMeshObject, ijkZone>(mesh),
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
  //  const faceList& faces=mesh_.faces();
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


  //if not in the Cartesian interface domain, do nothing.

  Info<<"Size of the interface domain of the master proc="<<cellIndices.size()<<endl;
  Pmin_=point::uniform(VGREAT);
  Pmax_=point::uniform(-VGREAT);
  dx_=0.0;dy_=0.0;dz_=0.0;


  if (cellIndices.size()>0)
    {

      ///// Old dx,dy, and dz calculation ////
      //  
      //  const label iCell=cellIndices[0];
      //  for (int iF=0;iF<6;iF++)
      // {	    
      //   label f1=cells[iCell][iF];
      //   label f2=cells[iCell].opposingFaceLabel(f1,faces);
      //   scalar dx=Cf[f1].x()-Cf[f2].x();
      //  scalar dy=Cf[f1].y()-Cf[f2].y();
      //   scalar dz=Cf[f1].z()-Cf[f2].z();
      //   if (dx>dx_)
      //	dx_=dx;
      //   if (dy>dy_)//
      //  dy_=dy;
      //   if (dz>dz_
      //	dz_=dz;		    
      //  }
      ///////////////////////////////////////
      const label iCell=cellIndices[0];
      const cell& c = mesh_.cells()[iCell];
      const boundBox& bb = c.localPoints(mesh_);
      dx_ = bb.span().x();
      dy_ = bb.span().y();
      dz_ = bb.span().z();
      ////////////////////////////////////

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
  

  ///////////////////// REFACTORING THE PARALLEL COMMUNICATION ///////////////////////
  if (cellIndices.size()>0)
    {
      const polyBoundaryMesh& boundary = mesh.boundaryMesh();
      HashSet<label> neighbourProcs;

      forAll(boundary, patchi)
	{
	  const polyPatch& patch = boundary[patchi];
	
	  if (isA<processorPolyPatch>(patch))
	    {
	      const processorPolyPatch& procPatch = refCast<const processorPolyPatch>(patch);
	      label nbrProcNo = procPatch.neighbProcNo();
	      neighbourProcs.insert(nbrProcNo);
	    }
	}

      // Output:
      Info << "This processor has neighbors: " << neighbourProcs.toc() << endl;

      labelList globalIds(cellIndices.size(),-10);

      
      forAll(cellIndices, i)
      {
	const label iCell = cellIndices[i];
	label ijk = ijk1(ijk3(iCell));
	globalIds[ijk] = globalNumbering_.toGlobal(iCell);
      }
    
      if (Pstream::parRun())
	{

	  ////// start from here ......

	  labelList myIds;
	  forAll(cellIndices, i)
	    {
	      const label iCell = cellIndices[i];
	      label ijk = ijk1(ijk3(iCell));
	      myIds.append(ijk);
	      myIds.append(globalNumbering_.toGlobal(iCell));
	    }
	
	  List<labelList> allMyIds(Pstream::nProcs());
	  allMyIds[Pstream::myProcNo()] = myIds;
	  Pstream::gatherList(allMyIds);
	  Pstream::broadcastList(allMyIds);
	
	  for (label procI = 0; procI < Pstream::nProcs(); ++procI)
	    {
	      if (procI != Pstream::myProcNo())
		{
		  const labelList& procIds = allMyIds[procI];
		  for (label j = 0; j < procIds.size(); j += 2)
		    {
		      globalIds[procIds[j]] = procIds[j+1];
		    }
		}
	    }
	}


    ///////////////////////////////////////////////////////////////////////////////////////




    
    // 1. Define the path to the file using the correct function.
    //    mesh_.time().casePath() ALWAYS points to the top-level case directory.
    fileName globalIdsPath = mesh_.time().casePath()/"constant"/"globalIds.dat";

    // 2. Check if the master file exists on disk.
    if (isFile(globalIdsPath))
    {
        // --- PATH 1: FILE EXISTS ---
        Info<< "Found pre-computed globalIds.dat. Reading from disk..." << endl;

        // Create an IOobject specifically for reading the master file.
        IOobject header
        (
            "globalIds.dat",
            mesh_.time().constant(), // The registry path is still the local one
            "constant",              // The instance is the top-level constant dir
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            true // Allow reading from outside the processor directory
        );
        // Perform the read operation directly into the member variable.
        header.read(globalIds_);
        Info<< "Successfully read " << globalIds_.size() << " entries." << endl;
    }
    else
    {
        // --- PATH 2: FILE DOES NOT EXIST ---
        // This part is unchanged. It generates the map on the fly
        // using the scalable gather/broadcast pattern and then writes
        // the file for the master processor.
        Info<< "Could not find globalIds.dat. Generating on the fly..." << endl;

        labelList globalIds(Nx_*Ny_*Nz_, -10);

        forAll(cellIndices, i)
        {
            const label iCell = cellIndices[i];
            label ijk = ijk1(ijk3(iCell));
            globalIds[ijk] = globalNumbering_.toGlobal(iCell);
        }

        if (Pstream::parRun())
        {
            labelList myIds;
            forAll(cellIndices, i)
            {
                const label iCell = cellIndices[i];
                label ijk = ijk1(ijk3(iCell));
                myIds.append(ijk);
                myIds.append(globalNumbering_.toGlobal(iCell));
            }

            List<labelList> allMyIds(Pstream::nProcs());
            allMyIds[Pstream::myProcNo()] = myIds;
            Pstream::gatherList(allMyIds);
            Pstream::broadcastList(allMyIds);

            for (label procI = 0; procI < Pstream::nProcs(); ++procI)
            {
                if (procI != Pstream::myProcNo())
                {
                    const labelList& procIds = allMyIds[procI];
                    for (label j = 0; j < procIds.size(); j += 2)
                    {
                        globalIds[procIds[j]] = procIds[j+1];
                    }
                }
            }
        }

        globalIds_ = globalIds;

        if (Pstream::master())
        {
            IOobject writeHeader
            (
                "globalIds.dat",
                mesh_.time().constant(),
                "constant",
                mesh_,
                IOobject::NO_READ,
                IOobject::WRITE_IF_PRESENT,
                true
            );
            Ostream& os = writeHeader.writeStream(IOstream::ASCII);
            os << globalIds_;
            os.flush();
            Info<< "Wrote newly generated globalIds.dat to disk for future runs." << endl;
        }
    }    
    

    Info<<"is x periodic: " <<isCyclic_.x()<<endl;
    Info<<"is y periodic: " <<isCyclic_.y()<<endl;
    Info<<"is z periodic: " <<isCyclic_.z()<<endl;

    //mark boundary cells
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
 


// ************************************************************************* //
