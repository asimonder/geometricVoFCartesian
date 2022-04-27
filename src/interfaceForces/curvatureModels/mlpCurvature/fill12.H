K_.correctBoundaryConditions();
Kf_=fvc::interpolate(K_);
const scalar d=2.;
neiMax=1;
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
	      point cc = C[iCell];
	      if (!interfaceCells[iCell])
		{
		  neiInterface[iCell]=true;
		}
	    }
	}
    }
}

alphaIJK.clear();
ijkMesh_.getZoneField(neiInterface,alphaIJK,alpha1_,stencilSize);
ijkMesh_.getZoneField(neiInterface,curvIJK,K_,stencilSize);
//ijkMesh_.getZoneField(neiInterface,faceCentreIJK,faceCentre,stencilSize);
const scalar distMax=sqrt(1.25)*dl_;
const surfaceVectorField& Cf = mesh.Cf();
const surfaceVectorField& Sf = mesh.Sf();
const surfaceScalarField& magSf = mesh.magSf();
const Foam::labelList& nei=alpha1_.mesh().faceNeighbour();
const Foam::labelList& own=alpha1_.mesh().faceOwner();
Info<<"nei="<<nei.size()<<endl;
Info<<"own="<<own.size()<<endl;
Info<<"Kf_.internal="<<Kf_.internalField().size()<<endl;  
Info<<"mesh.boundaryMesh="<<mesh.boundaryMesh().size()<<endl;
forAll(Kf_,iFace)
{
  label iCellO=own[iFace];
  label iCellN=nei[iFace];
  scalar K1=K_[own[iFace]];      
  scalar K2=K_[nei[iFace]];
  if (Kf_[iFace]!=0) // and K1*K2==0)
    {
      Kf_[iFace]=0;
      const point Cfp=Cf[iFace];
	      
      label iCell;
      if (interfaceCells[iCellO])
	iCell=iCellO;
      else
	iCell=iCellN;

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
		  if ((iP+i+1>Nx and ijkMesh_.symXOut()) or (iP+i<0 and ijkMesh_.symXIn())) 
		    il=iP-i;
		  else
		    il=(iP+i+Nx)%Nx;
		  if ((jP+j+1>Ny and ijkMesh_.symYOut()) or (jP+j<0 and ijkMesh_.symYIn())) 
		    jl=jP-j;
		  else
		    jl=(jP+j+Ny)%Ny;
		  if ((kP+k+1>Nz and ijkMesh_.symZOut()) or (kP+k<0 and ijkMesh_.symZIn())) 
		    kl=kP-k;
		  else
		    kl=(kP+k+Nz)%Nz;
		  label ijk=il+Nx*jl+Nx*Ny*kl;
		  label gblId=globalIds[ijk];
		  //vector faceCentre=faceCentreIJK[gblId];
		  //scalar dist = mag(Cfp-faceCentre);
		  if (alphaIJK[gblId]>interfaceTol_ && alphaIJK[gblId]<1-interfaceTol_)
		    {
		      scalar weight=1./Foam::pow(distCell,d);
		      sumW+= weight;
		      Kf_[iFace]+=weight*curvIJK[gblId];
		    }
		}
	    }
	}
      Kf_[iFace]/=sumW;
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
	  if (Kf_.boundaryFieldRef()[iPatch][iFace]!=0)
	    {
	      scalar Kff=Kf_.boundaryFieldRef()[iPatch][iFace];
	      Kf_.boundaryFieldRef()[iPatch][iFace]=0;
	      const point Cfp=Cf.boundaryField()[iPatch][iFace];
	      const vector Sfv=Sf.boundaryField()[iPatch][iFace]/magSf.boundaryField()[iPatch][iFace];
	      
	      point cc = C[iCell];
	      label iP,jP,kP;
	      if (!interfaceCells[iCell])
		{
		  cc=cc+2*(Cfp-cc);
		  //iP=round((cc.x()-Pmin.x())/dl_);
		  //jP=round((cc.y()-Pmin.y())/dl_);
		  //kP=round((cc.z()-Pmin.z())/dl_);
		  //printf("Not a interface cell in proc=%d. cc=(%f,%f,%f); i=(%d,%d,%d).\n",Pstream::myProcNo(),cc.x(),cc.y(),cc.z(),iP,jP,kP);
		  //cc=cc+2*(Cfp-cc);
		  //label iP2=round((cc.x()-Pmin.x())/dl_);
		  //label jP2=round((cc.y()-Pmin.y())/dl_);
		  //label kP2=round((cc.z()-Pmin.z())/dl_);
		  //printf("Not a interface cell in proc=%d. pType=%s; Kf=%f; Sf=(%.1f,%.1f,%.1f); i1=(%d,%d,%d); i2=(%d,%d,%d).\n",Pstream::myProcNo(),pType,Kff,Sfv.x(),Sfv.y(),Sfv.z(),iP,jP,kP,iP2,jP2,kP2);
		  //continue;
		}
	      iP=round((cc.x()-Pmin.x())/dl_);
	      jP=round((cc.y()-Pmin.y())/dl_);
	      kP=round((cc.z()-Pmin.z())/dl_);

	      if (iP+1>Nx or iP<0) 
		iP=(iP+Nx)%Nx;
	      if (jP+1>Ny or jP<0) 
		jP=(jP+Ny)%Ny;
	      if (kP+1>Nz or kP<0) 
		kP=(kP+Nz)%Nz;
	      
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
			  if ((iP+i+1>Nx and ijkMesh_.symXOut()) or (iP+i<0 and ijkMesh_.symXIn())) 
			    il=iP-i;
			  else
			    il=(iP+i+Nx)%Nx;
			  if ((jP+j+1>Ny and ijkMesh_.symYOut()) or (jP+j<0 and ijkMesh_.symYIn())) 
			    jl=jP-j;
			  else
			    jl=(jP+j+Ny)%Ny;
			  if ((kP+k+1>Nz and ijkMesh_.symZOut()) or (kP+k<0 and ijkMesh_.symZIn())) 
			    kl=kP-k;
			  else
			    kl=(kP+k+Nz)%Nz;
			  label ijk=il+Nx*jl+Nx*Ny*kl;
			  label gblId=globalIds[ijk];
			  //vector faceCentre=faceCentreIJK[gblId];
			  //scalar dist = mag(Cfp-faceCentre);
			  if (alphaIJK[gblId]>interfaceTol_ && alphaIJK[gblId]<1-interfaceTol_)
			    {
			      scalar weight=1./Foam::pow(distCell,d);
			      sumW+= weight;
			      Kf_.boundaryFieldRef()[iPatch][iFace]+=weight*curvIJK[gblId];
			    }
			}
		    }
		}
	      Kf_.boundaryFieldRef()[iPatch][iFace]/=sumW;
	    }
	}
    }
}  

