const scalar d=2.;
neiMax=1;
stencilSize=Vector<label>(neiMax,neiMax,neiMax);
K_.correctBoundaryConditions();
curvIJK.clear();
ijkMesh_.getZoneField(neiInterface,curvIJK,K_,stencilSize);
ijkMesh_.getZoneField(neiInterface,faceCentreIJK,faceCentre,stencilSize);

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

const scalar Lx=ijkMesh_.Lx()+dl_;
const scalar Ly=ijkMesh_.Ly()+dl_;
const scalar Lz=ijkMesh_.Lz()+dl_;
forAll(neiInterface,iCell)
{
  if ((alpha1_[iCell]<interfaceTol_ or alpha1_[iCell]>1-interfaceTol_) and interfaceCells[iCell])
    {
      const point cc = C[iCell];
      
      label iP=round((cc.x()-Pmin.x())/dl_);
      label jP=round((cc.y()-Pmin.y())/dl_);
      label kP=round((cc.z()-Pmin.z())/dl_);
      
      scalar sumW=0;
      label count=0;
      for (int i=-iMax;i<iMax+1;i++)
	{
	  for (int j=-jMax;j<jMax+1;j++)
	    {
	      for (int k=-kMax;k<kMax+1;k++)
		{
		  point pos(dl_*i,dl_*j,dl_*k);
		  scalar distCell = mag(pos);
		  if (distCell>=1.99*dl_ or distCell<0.01*dl_)
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

		  label ijk=il+Nx*jl+Nx*Ny*kl;
		  label gblId=globalIds[ijk];
		  if (curvIJK[gblId]!=0)
		    {
		      vector fC=faceCentreIJK[gblId];
		      vector dist = cc-fC;
		      if (dist.x()>2*dl_)
			dist.x()-=Lx;
		      else if (dist.x()<-2*dl_)
			dist.x()+=Lx;
		      if (dist.y()>2*dl_)
			dist.y()-=Ly;
		      else if (dist.y()<-2*dl_)
			dist.y()+=Ly;
		      if (dist.z()>2*dl_)
			dist.z()-=Lz;
		      else if (dist.z()<-2*dl_)
			dist.z()+=Lz;
		      scalar distI= Foam::mag(dist);
		      if (distI>2*dl_)
			FatalErrorInFunction
			  << "Extrapolation to neightbouring cell has failed! Distance is too large. d="<<distI
			  << abort(FatalError);
		      scalar weight=1./Foam::pow(distI,d);
		      sumW+= weight;
		      count+=1;
		      K_[iCell]+=weight*curvIJK[gblId];
		    }
		}
	    }
	}
      if (sumW>0)
	{
	  K_[iCell]/=sumW;
	}
    }
}
