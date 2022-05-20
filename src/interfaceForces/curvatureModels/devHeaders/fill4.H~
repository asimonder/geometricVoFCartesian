neiMax=2;
stencilSize=Vector<label>(neiMax,neiMax,neiMax);

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

alphaIJK.clear();
curvIJK.clear();
ijkMesh_.getZoneField(neiInterface,alphaIJK,alpha1_,stencilSize);
ijkMesh_.getZoneField(neiInterface,faceCentreIJK,faceCentre,stencilSize);
ijkMesh_.getZoneField(neiInterface,curvIJK,K_,stencilSize);

forAll(neiInterface,celli)
{
  if (neiInterface[celli])
    {
      const point cc = C[celli];
      
      label iP=round((cc.x()-Pmin.x())/dl_);
      label jP=round((cc.y()-Pmin.y())/dl_);
      label kP=round((cc.z()-Pmin.z())/dl_);
      
      label gblIdMin=-1;
      scalar distMinCell=GREAT;
      scalar distMinPlane=GREAT;
      for (int i=-iMax;i<iMax+1;i++)
	{
	  for (int j=-jMax;j<jMax+1;j++)
	    {
	      for (int k=-kMax;k<kMax+1;k++)
		{
		  label sum=mag(i)+mag(j)+mag(k);
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
		  vector faceCentre=faceCentreIJK[gblId];
		  scalar distPlane = mag(cc-faceCentre);
		  vector pos(dl_*i,dl_*j,dl_*k);
		  scalar distCell = mag(pos);
		  if (mag(faceCentre)>0 && alphaIJK[gblId]>interfaceTol_ && alphaIJK[gblId]<1-interfaceTol_)
		    {
		      if (distMinCell-distCell>0.0001*dl_) // and distPlane < distMinPlane)
			{
			  distMinCell = distCell;
			  distMinPlane = distPlane;
			  gblIdMin = gblId;
			}
		      else if (mag(distCell-distMinCell)<0.0001*dl_  and distPlane < distMinPlane)
			{
			  distMinCell = distCell;
			  distMinPlane = distPlane;
			  gblIdMin = gblId;
			}
		    }
		}
	    }
	}
      if (gblIdMin!=-1)
	{
	  K_[celli]=curvIJK[gblIdMin];
	}
    }
}
