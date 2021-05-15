/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2021 Asim Onder
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

  fileName dictPath=mesh_.time().path()/"system/ijkDict";
  if (Pstream::parRun())
    dictPath=mesh_.time().path()/"../system/ijkDict";

  IFstream is(dictPath);
  dictionary ijkDict(is);

  word cellSetName=ijkDict.get<word>("cellSet");
  bool isDomain=true;
  if (cellSetName!="domain")
    isDomain=false;

  const  surfaceVectorField&  Cf = mesh_.Cf();

  cellList cells=mesh_.cells();
  faceList faces=mesh_.faces();
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

  Info<<"Size of the interface domain="<<cellIndices.size()<<endl;
  dx_=0.0;dy_=0.0;dz_=0.0;
  Pmin_=point::uniform(VGREAT);
  Pmax_=point::uniform(-VGREAT);
  if (cellIndices.size()>0)
    {
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
      //for(int celli=0;celli<100<celli++)
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
    
    labelList globalIds(Nx_*Ny_*Nz_,-10);
  
    forAll(cellIndices,celli)
      {
	const label iCell=cellIndices[celli];
	//Vector<label> ijk3=ijk3(iCell);
	label ijk=ijk1(ijk3(iCell));
	if (Pstream::parRun())
	  globalIds[ijk]=globalNumbering_.toGlobal(iCell);
	else
	  globalIds[ijk]=iCell;
      }
    
    if (Pstream::parRun())
      {
	Foam::reduce(globalIds, maxOp<labelList>());
      }

    globalIds_=globalIds;
    
    //isCyclic?
    forAll(faces,facei)
      {
	if (!mesh_.isInternalFace(facei))
	  {
	    if ((C_[mesh_.faceOwner()[facei]]<Pmin_ && C_[mesh_.faceNeighbour()[facei]]<Pmin_) ||
		(C_[mesh_.faceOwner()[facei]]>Pmax_ && C_[mesh_.faceNeighbour()[facei]]>Pmax_))
	      {
		continue;
	      }
	    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();
	    	    
	    // Boundary face. Find out which face of which patch
	    const label patchi = pbm.patchID()[facei - mesh_.nInternalFaces()];
    
	    // Handle empty patches
	    const polyPatch& pp = pbm[patchi];
	    //printf("Proc=%d, patch type=%s \n",Pstream::myProcNo(),mesh_.boundary()[patchi].type());
	    if (isA<cyclicPolyPatch>(pp) || mesh_.boundary()[patchi].type()=="processorCyclic")
	      {
		//printf("Proc=%d, cyclic patch found...\n",Pstream::myProcNo());
		vector transDir=mesh_.faceAreas()[facei];
		if (transDir.x()!=0 and transDir.y()==0 and transDir.z()==0)
		  isCyclic_.x()=true;
		else if (transDir.y()!=0 and transDir.x()==0 and transDir.z()==0)
		  isCyclic_.y()=true;
		else if (transDir.z()!=0 and transDir.x()==0 and transDir.y()==0)
		  isCyclic_.z()=true;
		else
		  FatalErrorInFunction
		    << "Cyclic BC found but not aligned with one of the Cartesian coordinates"
		    << abort(FatalError);
	      }
	  }
      }

    //printf("Proc=%d, isCyclic.x()=%d, isCyclic.y()=%d, isCyclic.z()=%d\n",Pstream::myProcNo(),isCyclic_.x(),isCyclic_.y(),isCyclic_.z());
    
    if (Pstream::parRun())
      {
	Foam::reduce(isCyclic_.x(),orOp<bool>());
	Foam::reduce(isCyclic_.y(),orOp<bool>());
	Foam::reduce(isCyclic_.z(),orOp<bool>());
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
      if ((!isCyclic_.x() and !isEmpty_.x() and (ijk.x()+nOffsetCells>=Nx_ or ijk.x()-nOffsetCells<=0)) or
	  (!isCyclic_.y() and !isEmpty_.y() and (ijk.y()+nOffsetCells>=Ny_ or ijk.y()-nOffsetCells<=0)) or
	  (!isCyclic_.z() and !isEmpty_.z() and (ijk.z()+nOffsetCells>=Nz_ or ijk.z()-nOffsetCells<=0)))
	atBoundary[celli]=true;
    }

  //return atBoundary;      

}


Foam::Vector<Foam::label> Foam::ijkZone::ijk3(label celli)
{  
  label i=round((C_[celli].x()-Pmin_.x())/dx_);
  label j=round((C_[celli].y()-Pmin_.y())/dy_);
  label k=round((C_[celli].z()-Pmin_.z())/dz_);

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
