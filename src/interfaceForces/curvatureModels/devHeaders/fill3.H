Info<<"Neighbouring cells of the interface are left zero."<<endl;
forAll(interfaceCells,iCell)
{
  if(interfaceCells[iCell])
    {
      neiInterface[iCell]=true;
      if (alpha1_[iCell]<interfaceTol_ or alpha1_[iCell]>1-interfaceTol_) //mag(faceCentre)>
	{
	    K_[iCell]=0.0;
	}
    }
}

if (extendInterface_)
  {
#include "extendInterface3.H"
  }

forAll(mesh.boundaryMesh(),iPatch)
{
  word pType=Kf_.boundaryField()[iPatch].type();
  if(pType!="empty")
    {	  
      const polyPatch& cPatch = mesh.boundaryMesh()[iPatch];
      const labelUList& faceCells = cPatch.faceCells();
      forAll(faceCells,iFace)
	{
	  label iCell=faceCells[iFace];
	  if (Kf_.boundaryFieldRef()[iPatch][iFace]!=0)
	    {
	      neiInterface[iCell]=true;
	    }
	}
    }
}

K_.correctBoundaryConditions();
curvIJK.clear();
ijkMesh_.getZoneField(neiInterface,curvIJK,K_,stencilSize);

const surfaceVectorField& Cf = mesh.Cf();
Kf_=fvc::interpolate(K_);
const Foam::labelList& nei=alpha1_.mesh().faceNeighbour();
const Foam::labelList& own=alpha1_.mesh().faceOwner();
forAll(Kf_,iFace)
{
  if (Kf_[iFace]!=0)
    {
      scalar K1=K_[own[iFace]];      
      scalar K2=K_[nei[iFace]];
      if (K1*K2==0)
	{
	  if (mag(K1)>mag(K2))
	    Kf_[iFace]=K1;
	  else
	    Kf_[iFace]=K2;
	}
    }
}

forAll(mesh.boundaryMesh(),iPatch)
{
  word pType=Kf_.boundaryField()[iPatch].type();
  if(pType!="empty")
    {	  
      const polyPatch& cPatch = mesh.boundaryMesh()[iPatch];
      label iFaceStart = cPatch.start();
      const labelUList& faceCells = cPatch.faceCells();
      forAll(faceCells,iFace)
	{
	  label iCell=faceCells[iFace];
	  if(neiInterface[iCell])
	    {
	      const point cc = C[iCell];
	      label iP=round((cc.x()-Pmin.x())/dl_);
	      label jP=round((cc.y()-Pmin.y())/dl_);
	      label kP=round((cc.z()-Pmin.z())/dl_);
	      label ijk=iP+Nx*jP+Nx*Ny*kP;
	      label gblId=globalIds[ijk];
	      scalar Kp=curvIJK[gblId];
	      const point cfp=Cf.boundaryField()[iPatch][iFace];
	      const point R = cfp-cc;
	      label iN,jN,kN;
		  
	      iN=round(2.0*R.x()/dl_);
	      jN=round(2.0*R.y()/dl_);
	      kN=round(2.0*R.z()/dl_);
	      //Info<<"R="<<R<<";r=("<<iN<<","<<jN<<","<<kN<<")"<<endl;				  
	      if(pType=="symmetry" or pType=="wall")
		continue; //iN=iP-iN;jN=jP-jN;kN=kP-kN;
	      else //if (pType=="cyclic")
		{
		  iP=(iP+iN+Nx)%Nx;
		  jP=(jP+jN+Ny)%Ny;
		  kP=(kP+kN+Nz)%Nz;
		}
	      ijk=iP+Nx*jP+Nx*Ny*kP;
	      gblId=globalIds[ijk];
	      scalar Kn=curvIJK[gblId];
	      if (Kn*Kp==0)
		{
		  if (mag(Kn)>mag(Kp))
		    Kf_.boundaryFieldRef()[iPatch][iFace]=Kn;
		  else
		    Kf_.boundaryFieldRef()[iPatch][iFace]=Kp;
		  //if(Kf_.boundaryField()[iPatch][iFace]!=0)
		  //	printf("Patch: %s; Kf(%f,%f,%f)=%f; iN=(%d,%d,%d); in proc=%d\n",pType,cfp.x(),cfp.y(),cfp.z()
		  //       ,Kf_.boundaryField()[iPatch][iFace],iN,jN,kN,Pstream::myProcNo());
		  //Info<<"Patch: "<<pType<<" ;Kf_.boundaryField()[iPatch][iFace]="<<Kf_.boundaryField()[iPatch][iFace]<<endl;
		}
	    }
	}
    }
}
