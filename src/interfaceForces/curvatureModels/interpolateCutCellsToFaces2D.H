const volVectorField& faceCentre = surf.centre();
point Pmin=ijkMesh_.Pmin();
//point Pmax=ijkMesh_.Pmax();
//const scalar dx=ijkMesh_.dx();
//const scalar dy=ijkMesh_.dy();
//const scalar dz=ijkMesh_.dz();
const label Nx=ijkMesh_.Nx();
const label Ny=ijkMesh_.Ny();
const label Nz=ijkMesh_.Nz();
Map<vector> faceCentreIJK;
Map<scalar> curvIJK;

K_.correctBoundaryConditions();
Kf_=fvc::interpolate(K_);

const scalar Lx=ijkMesh_.Lx()+dl_;
const scalar Ly=ijkMesh_.Ly()+dl_;
const scalar Lz=ijkMesh_.Lz()+dl_;

label neiMax=1;
stencilSize=Vector<label>(1,1,1);
label iMax,jMax,kMax;
label il,jl,kl;
boolList neiInterface(mesh.nCells(),false);

const scalar d=2.;
if (ijkMesh_.isEmpty().x())
  {
    iMax=0;jMax=neiMax;kMax=neiMax;
  }
 else if (ijkMesh_.isEmpty().y())
   {
     iMax=neiMax;jMax=0;kMax=neiMax;
   }
 else if (ijkMesh_.isEmpty().z())
   {
     iMax=neiMax;jMax=neiMax;kMax=0;
   }
 else
   {
     iMax=neiMax;jMax=neiMax;kMax=neiMax;
   }

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

/*if (extendInterface_)
  {
#include "extendInterface3.H"
  }
*/

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
	      if (!interfaceCells[iCell])
		{
		  neiInterface[iCell]=true;
		}
	    }
	}
    }
}


alphaIJK.clear();
//ijkMesh_.getZoneField(neiInterface,alphaIJK,alpha1_,stencilSize);
ijkMesh_.getZoneField(neiInterface,curvIJK,K_,stencilSize);
Info<<"Done: ijkMesh_.getZoneField(neiInterface,curvIJK,K_,stencilSize);"<<endl;
ijkMesh_.getZoneField(neiInterface,faceCentreIJK,faceCentre,stencilSize);
Info<<"Done: ijkMesh_.getZoneField(neiInterface,faceCentreIJK,faceCentre,stencilSize);"<<endl;
const scalar distMax=1.01*sqrt(1.25)*dl_;
const scalar distInterfaceMax=4.0*dl_; 
const surfaceVectorField& Cf = mesh.Cf();
const Foam::labelList& nei=alpha1_.mesh().faceNeighbour();
const Foam::labelList& own=alpha1_.mesh().faceOwner();
forAll(Kf_,iFace)
{
  label iCellO=own[iFace];
  scalar K1=K_[own[iFace]];      
  scalar K2=K_[nei[iFace]];
  if (Kf_[iFace]!=0 and K1*K2==0)
    {
      Kf_[iFace]=0;
      const point Cfp=Cf[iFace];
	      
      label iCell=iCellO;
      point cc = C[iCell];
      label iP=round((cc.x()-Pmin.x())/dl_);
      label jP=round((cc.y()-Pmin.y())/dl_);
      label kP=round((cc.z()-Pmin.z())/dl_);
      
      scalar sumW=0;
      for (int i=-iMax;i<iMax+1;i++)
	{
	  for (int j=-jMax;j<jMax+1;j++)
	    {
	      for (int k=-kMax;k<kMax+1;k++)
		{
		  point pos(dl_*i,dl_*j,dl_*k);
		  pos+=cc;
		  scalar distCell = mag(pos-Cfp);
		  if (distCell>distMax)
		    continue;
		  if (iP+i+1>Nx and ijkMesh_.symXOut())
		    il=iP+1-i;
		  else if (iP+i<0 and ijkMesh_.symXIn())
		    il=iP-1-i;
		  else
		    il=(iP+i+Nx)%Nx;
		  if (jP+j+1>Ny and ijkMesh_.symYOut()) 
		    jl=jP+1-j;
		  else if (jP+j<0 and ijkMesh_.symYIn())
		    jl=jP-1-j;
		  else
		    jl=(jP+j+Ny)%Ny;
		  if (kP+k+1>Nz and ijkMesh_.symZOut()) 
		    kl=kP+1-k;
		  else if (kP-k<0 and ijkMesh_.symZIn())
		    kl=kP-1-k;
		  else
		    kl=(kP+k+Nz)%Nz;
		  label gblId=ijkMesh_.lookupGlobalId(il,jl,kl);
		  
		  if (curvIJK[gblId]!=0)
		    {
		      vector fC=faceCentreIJK[gblId];
		      vector dist = fC-Cfp;
			      
		      if (dist.x()>distInterfaceMax)
			dist.x()-=Lx;
		      else if (dist.x()<-distInterfaceMax)
			dist.x()+=Lx;
		      if (dist.y()>distInterfaceMax)
			dist.y()-=Ly;
		      else if (dist.y()<-distInterfaceMax)
			dist.y()+=Ly;
		      if (dist.z()>distInterfaceMax)
			dist.z()-=Lz;
		      else if (dist.z()<-distInterfaceMax)
			dist.z()+=Lz;
		      scalar distI= Foam::mag(dist);
		      if (distI>3*dl_)
			{
			  FatalErrorInFunction
			    << "Extrapolation to neightbouring cell has failed! Distance is too large. i="<<i<<"; j="<<j<<"; k="<<k<<"; Cfp="<<Cfp<<"; fC="<<fC<<"; distCell/dl="<<distCell/dl_<<"; distI/dl="<<distI/dl_<<"; dist="<<dist/dl_<< abort(FatalError);
			}
		      scalar weight=1./Foam::pow(distI,d);
		      sumW+= weight;
		      Kf_[iFace]+=weight*curvIJK[gblId];
		    }
		}
	    }
	}
      if (sumW>0)
	{
	  Kf_[iFace]/=sumW;
	  //printf("Face (%f,%f) in proc=%d. Kf_=%.2f. count=%d. sum of weights=%.2f\n",Cfp.x(),Cfp.z(),Pstream::myProcNo(),Kf_[iFace],count,sumW);
	}

    }
}

forAll(mesh.boundaryMesh(),iPatch)
{
  word pType=Kf_.boundaryField()[iPatch].type();
  if(pType!="empty")
    {	  
      const polyPatch& cPatch = mesh.boundaryMesh()[iPatch];
      //label iFaceStart = cPatch.start();
      const labelUList& faceCells = cPatch.faceCells();
      forAll(faceCells,iFace)
	{
	  if (Kf_.boundaryFieldRef()[iPatch][iFace]!=0)
	    {
	      Kf_.boundaryFieldRef()[iPatch][iFace]=0;
	      const point Cfp=Cf.boundaryField()[iPatch][iFace];

	      label iCell=faceCells[iFace];
	      point cc = C[iCell];
	      label iP,jP,kP;
	      iP=round((cc.x()-Pmin.x())/dl_);
	      jP=round((cc.y()-Pmin.y())/dl_);
	      kP=round((cc.z()-Pmin.z())/dl_);
	      
	      scalar sumW=0;
	      for (int i=-iMax;i<iMax+1;i++)
		{
		  for (int j=-jMax;j<jMax+1;j++)
		    {
		      for (int k=-kMax;k<kMax+1;k++)
			{
			  point pos(dl_*i,dl_*j,dl_*k);
			  pos+=cc;
			  scalar distCell = mag(pos-Cfp);
			  if (distCell>distMax)
			    continue;			  
			  if (iP+i+1>Nx and ijkMesh_.symXOut())
			    il=iP+1-i;
			  else if (iP+i<0 and ijkMesh_.symXIn())
			    il=iP-1-i;
			  else
			    il=(iP+i+Nx)%Nx;
			  if (jP+j+1>Ny and ijkMesh_.symYOut()) 
			    jl=jP+1-j;
			  else if (jP+j<0 and ijkMesh_.symYIn())
			    jl=jP-1-j;
			  else
			    jl=(jP+j+Ny)%Ny;
			  if (kP+k+1>Nz and ijkMesh_.symZOut()) 
			    kl=kP+1-k;
			  else if (kP-k<0 and ijkMesh_.symZIn())
			    kl=kP-1-k;
			  else
			    kl=(kP+k+Nz)%Nz;
			  label gblId=ijkMesh_.lookupGlobalId(il,jl,kl);
			  //if (alphaIJK[gblId]>interfaceTol_ && alphaIJK[gblId]<1-interfaceTol_ && curvIJK[gblId]!=0)
			  if (curvIJK[gblId]!=0)
			    {
			      vector fC=faceCentreIJK[gblId];
			      vector dist = fC-Cfp;
			      
			      if (dist.x()>distInterfaceMax)
				dist.x()-=Lx;
			      else if (dist.x()<-distInterfaceMax)
				dist.x()+=Lx;
			      if (dist.y()>distInterfaceMax)
				dist.y()-=Ly;
			      else if (dist.y()<-distInterfaceMax)
				dist.y()+=Ly;
			      if (dist.z()>distInterfaceMax)
				dist.z()-=Lz;
			      else if (dist.z()<-distInterfaceMax)
				dist.z()+=Lz;
			      scalar distI= Foam::mag(dist);
			      if (distI>3*dl_)
				FatalErrorInFunction
				  << "Extrapolation to neightbouring cell has failed! Distance is too large. fC="<<fC<<"; distCell="<<distCell<<"; distI="<<distI<<"; dist="<<dist<< abort(FatalError);
			      scalar weight=1./Foam::pow(distI,d);
			      sumW+= weight;
			      Kf_.boundaryFieldRef()[iPatch][iFace]+=weight*curvIJK[gblId];
			    }
			}
		    }
		}
	      if (sumW>0)
		{
		  Kf_.boundaryFieldRef()[iPatch][iFace]/=sumW;
		  //printf("Boundary face (%f,%f) in proc=%d. Kf_=%.2f. pType=%s. count=%d. sum of weights=%.2f\n",Cfp.x(),Cfp.z(),Pstream::myProcNo(),Kf_.boundaryFieldRef()[iPatch][iFace],pType,count, sumW);
		}
	    }
	}
    }
}  

 
