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

#include "uniformStencil.H"
#include <cassert>
// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

/*namespace Foam
{
  defineTypeNameAndDebug(uniformStencil, 0);
  }*/


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Protected constructor for subclasses (no ijkZone needed)
Foam::uniformStencil::uniformStencil
(
    SubclassTag,
    const fvMesh& mesh,
    int iMax,
    int jMax,
    int kMax,
    const Vector<bool>& isEmpty,
    scalar dx
)
:
    mesh_(mesh),
    ijkMeshPtr_(nullptr),
    N_(0),
    is2D_(false),
    iMax_(iMax),
    jMax_(jMax),
    kMax_(kMax),
    isEmpty_(isEmpty),
    dx_(dx)
{
    if (isEmpty_.x()) { iMax_ = 0; is2D_ = true; }
    if (isEmpty_.y()) { jMax_ = 0; is2D_ = true; }
    if (isEmpty_.z()) { kMax_ = 0; is2D_ = true; }

    NxS_ = 2 * iMax_ + 1;
    NyS_ = 2 * jMax_ + 1;
    NzS_ = 2 * kMax_ + 1;
    N_ = NxS_ * NyS_ * NzS_;

    A_.setSize(N_, 0.0);
}


Foam::uniformStencil::uniformStencil(const Foam::fvMesh& mesh,
				     const ijkZone& ijkMesh,
				     int iMax,
				     int jMax,
				     int kMax):
  mesh_(mesh),
  ijkMeshPtr_(&ijkMesh),
  is2D_(false),
  iMax_(iMax),
  jMax_(jMax),
  kMax_(kMax),
  isEmpty_(ijkMesh.isEmpty()),
  dx_(ijkMesh.dx())
{  

  if (isEmpty_.x())
    {
      iMax_=0;is2D_=true;
    }
  if (isEmpty_.y())
    {
      jMax_=0;is2D_=true;
    }
  if (isEmpty_.z())
    {
      kMax_=0;is2D_=true;
    }

  NxS_=2*iMax_+1;
  NyS_=2*jMax_+1;
  NzS_=2*kMax_+1;

  N_=NxS_*NyS_*NzS_;
  Info<<"uniformStencil: iMax_="<<iMax_<<endl;
  Info<<"uniformStencil: jMax_="<<jMax_<<endl;
  Info<<"uniformStencil: kMax_="<<kMax_<<endl;
  Info<<"uniformStencil: N_="<<N_<<endl;  
  A_.setSize(N_,0.0);  

}

Foam::uniformStencil::uniformStencil(const Foam::fvMesh& mesh,
				     const ijkZone& ijkMesh,
				     int nMax):
  mesh_(mesh),
  ijkMeshPtr_(&ijkMesh),
  is2D_(false),
  iMax_(nMax),
  jMax_(nMax),
  kMax_(nMax),
  NxS_(2*nMax+1),
  NyS_(2*nMax+1),
  NzS_(2*nMax+1),
  isEmpty_(ijkMesh.isEmpty()),
  dx_(ijkMesh.dx())
{  
  if (isEmpty_.x())
    {
      NxS_=1;is2D_=true;
    }
  if (isEmpty_.y())
    {
      NyS_=1;is2D_=true;
    }
  if (isEmpty_.z())
    {
      NzS_=1;is2D_=true;
    }

  Info<<"uniformStencil: evenly distributed construction..."<<endl;
  Info<<"uniformStencil even: iMax="<<iMax_<<endl;
  Info<<"uniformStencil even: jMax="<<jMax_<<endl;
  Info<<"uniformStencil even: kMax="<<kMax_<<endl;
  Info<<"uniformStencil even: is2D="<<is2D_<<endl;

  N_=NxS_*NyS_*NzS_;
  Info<<"uniformStencil even: N_="<<N_<<endl;
  A_.setSize(N_,0.0);  
}


void Foam::uniformStencil::setStencil(const Map<scalar>& phiIJK, const Vector<label>& ijk)
{
  const ijkZone& ijkMesh = *ijkMeshPtr_;
  const label Nx=ijkMesh.Nx();
  const label Ny=ijkMesh.Ny();
  const label Nz=ijkMesh.Nz();
  
  const label& iP=ijk.x();
  const label& jP=ijk.y();
  const label& kP=ijk.z();

  int il,jl,kl=0;
      
  if (isEmpty_.x())
    {
      il=0;
      for (int j=-jMax_;j<jMax_+1;j++)
	{
	  if (jP+j+1>Ny and ijkMesh.symYOut()) 
	    jl=jP+1-j;
	  else if (jP+j<0 and ijkMesh.symYIn())
	    jl=jP-1-j;
	  else
	    jl=(jP+j+Ny)%Ny;
	  for (int k=-kMax_;k<kMax_+1;k++)
	    {
	      if (kP+k+1>Nz and ijkMesh.symZOut()) 
		kl=kP+1-k;
	      else if (kP-k<0 and ijkMesh.symZIn())
		kl=kP-1-k;
	      else
		kl=(kP+k+Nz)%Nz;
	      label gblIdx=ijkMesh.lookupGlobalId(il,jl,kl);
	      A_[a2(j,k)]=phiIJK[gblIdx];	
	    }
	}
    }
  else if (isEmpty_.y())
    {
      jl=0;
      for (int i=-iMax_;i<iMax_+1;i++)
	{
	  if (iP+i+1>Nx and ijkMesh.symXOut())
	    il=iP+1-i;
	  else if (iP+i<0 and ijkMesh.symXIn())
	    il=iP-1-i;
	  else
	    il=(iP+i+Nx)%Nx;
	  for (int k=-kMax_;k<kMax_+1;k++)
	    {
	      if (kP+k+1>Nz and ijkMesh.symZOut()) 
		kl=kP+1-k;
	      else if (kP-k<0 and ijkMesh.symZIn())
		kl=kP-1-k;
	      else
		kl=(kP+k+Nz)%Nz;
	      label gblIdx=ijkMesh.lookupGlobalId(il,jl,kl);
	      A_[a2(i,k)]=phiIJK[gblIdx];
	    }
	}
    }
  else if (isEmpty_.z())
    {
      kl=0;
      for (int i=-iMax_;i<iMax_+1;i++)
	{
	  if (iP+i+1>Nx and ijkMesh.symXOut())
	    il=iP+1-i;
	  else if (iP+i<0 and ijkMesh.symXIn())
	    il=iP-1-i;
	  else
	    il=(iP+i+Nx)%Nx;
	  for (int j=-jMax_;j<jMax_+1;j++)
	    {
	      if (jP+j+1>Ny and ijkMesh.symYOut()) 
		jl=jP+1-j;
	      else if (jP+j<0 and ijkMesh.symYIn())
		jl=jP-1-j;
	      else
		jl=(jP+j+Ny)%Ny;
	      label gblIdx=ijkMesh.lookupGlobalId(il,jl,kl);
	      A_[a2(i,j)]=phiIJK[gblIdx];	
	    }
	}
    }
  else
    {
      for (int i=-iMax_;i<iMax_+1;i++)
	{
	  if (iP+i+1>Nx and ijkMesh.symXOut())
	    il=iP+1-i;
	  else if (iP+i<0 and ijkMesh.symXIn())
	    il=iP-1-i;
	  else
	    il=(iP+i+Nx)%Nx;
	  for (int j=-jMax_;j<jMax_+1;j++)
	    {
	      if (jP+j+1>Ny and ijkMesh.symYOut()) 
		jl=jP+1-j;
	      else if (jP+j<0 and ijkMesh.symYIn())
		jl=jP-1-j;
	      else
		jl=(jP+j+Ny)%Ny;
	      for (int k=-kMax_;k<kMax_+1;k++)
		{
		  if (kP+k+1>Nz and ijkMesh.symZOut()) 
		   kl=kP+1-k;
		 else if (kP-k<0 and ijkMesh.symZIn())
		   kl=kP-1-k;
		 else
		   kl=(kP+k+Nz)%Nz;
		  label gblIdx=ijkMesh.lookupGlobalId(il,jl,kl);
		  A_[a3(i,j,k)]=phiIJK[gblIdx];	
		}
	    }
	}
    }
}


Foam::vector Foam::uniformStencil::calcYoungNormal()
{
  Foam::vector m(0.0,0.0,0.0);

  const List<scalar>& A=A_;
  
  if (is2D_)
    {
      scalar mN=1./8./dx_*(A[a2(1,1)]-A[a2(-1,1)]+2.0*A[a2(1,0)]-2.0*A[a2(-1,0)]+A[a2(1,-1)]-A[a2(-1,-1)]);
      scalar mT=1./8./dx_*(A[a2(1,1)]-A[a2(1,-1)]+2.0*A[a2(0,1)]-2.0*A[a2(0,-1)]+A[a2(-1,1)]-A[a2(-1,-1)]);
      if (isEmpty_.x())
	m=vector(0,mN,mT);
      else if (isEmpty_.y())
	m=vector(mN,0,mT);
      else
	m=vector(mN,mT,0);
    }
  else
    {
      m.x()=1./32./dx_*(A[a3(1,-1,-1)]-A[a3(-1,-1,-1)]+2*A[a3(1,0,-1)]-2*A[a3(-1,0,-1)]+A[a3(1,1,-1)]-A[a3(-1,1,-1)]+
		      2*A[a3(1,-1,0)]-2*A[a3(-1,-1,0)]+4*A[a3(1,0,0)]-4*A[a3(-1,0,0)]+2*A[a3(1,1,0)]-2*A[a3(-1,1,0)]+
		      A[a3(1,-1,1)]-A[a3(-1,-1,1)]+2*A[a3(1,0,1)]-2*A[a3(-1,0,1)]+A[a3(1,1,1)]-A[a3(-1,1,1)]);

      m.y()=1./32./dx_*(A[a3(-1,1,-1)]-A[a3(-1,-1,-1)]+2*A[a3(0,1,-1)]-2*A[a3(0,-1,-1)]+A[a3(1,1,-1)]-A[a3(1,-1,-1)]+
		      2*A[a3(-1,1,0)]-2*A[a3(-1,-1,0)]+4*A[a3(0,1,0)]-4*A[a3(0,-1,0)]+2*A[a3(1,1,0)]-2*A[a3(1,-1,0)]+
		      A[a3(-1,1,1)]-A[a3(-1,-1,1)]+2*A[a3(0,1,1)]-2*A[a3(0,-1,1)]+A[a3(1,1,1)]-A[a3(1,-1,1)]);

      m.z()=1./32./dx_*(A[a3(-1,-1,1)]-A[a3(-1,-1,-1)]+2*A[a3(-1,0,1)]-2*A[a3(-1,0,-1)]+A[a3(-1,1,1)]-A[a3(-1,1,-1)]+
		      2*A[a3(0,-1,1)]-2*A[a3(0,-1,-1)]+4*A[a3(0,0,1)]-4*A[a3(0,0,-1)]+2*A[a3(0,1,1)]-2*A[a3(0,1,-1)]+
		      A[a3(1,-1,1)]-A[a3(1,-1,-1)]+2*A[a3(1,0,1)]-2*A[a3(1,0,-1)]+A[a3(1,1,1)]-A[a3(1,1,-1)]);
    }

  return m;
}

 
Foam::vector Foam::uniformStencil::calcCCDNormal() 
{

  const List<scalar>& A=A_;

  if (is2D_)
    {
      Foam::vector m(0.0,0.0,0.0);
      scalar mxc=0;
      scalar myc=0;
      scalar mxb=0;
      scalar myb=0;
      scalar mxf=0;
      scalar myf=0;
      scalar mx=0;
      scalar my=0;
      scalar sgn_mx=0;
      scalar sgn_my=0;

      for (int i=-1;i<2;i++)
	{
	  sgn_mx-=A[a2(1,i)]-A[a2(-1,i)];
	  sgn_my-=A[a2(i,1)]-A[a2(i,-1)];
	}

      if (sgn_mx!=0)
	sgn_mx/=fabs(sgn_mx);
      if (sgn_my!=0)
	sgn_my/=fabs(sgn_my);

      for (int i=-1;i<2;i++)
	{
	  mxc-=(A[a2(1,i)]-A[a2(-1,i)])/2.0;
	  myc-=(A[a2(i,1)]-A[a2(i,-1)])/2.0;
	  mxf-=(A[a2(1,i)]-A[a2(0,i)]);
	  myf-=(A[a2(i,1)]-A[a2(i,0)]);       
	  mxb-=(A[a2(0,i)]-A[a2(-1,i)]);
	  myb-=(A[a2(i,0)]-A[a2(i,-1)]);       
	}      

      mx=mxc;
      my=myc;
      
      if (fabs(my)<fabs(mx))
	{
	  mx=sgn_mx;
	}
      else
	{
	  my=sgn_my;
	}

      if (isEmpty_.x())
	m=vector(0,-mx,-my);
      else if (isEmpty_.y())
	m=vector(-mx,0,-my);
      else
	m=vector(-mx,-my,0);

      return m;
    }
  else
    {
      scalar sgn_mx=0;
      scalar sgn_my=0;
      scalar sgn_mz=0;
      scalar sgn_mxb=0;
      scalar sgn_myb=0;
      scalar sgn_mzb=0;

      for (int i=-1;i<2;i++)
	{
	  for (int j=-1;j<2;j++)
	    {
	      sgn_mx-=A[a3(1,i,j)]-A[a3(-1,i,j)];
	      sgn_my-=A[a3(i,1,j)]-A[a3(i,-1,j)];
	      sgn_mz-=A[a3(i,j,1)]-A[a3(i,j,-1)];
	    }
	}

      if (sgn_mx!=0)
	sgn_mx/=fabs(sgn_mx);
      if (sgn_my!=0)
	sgn_my/=fabs(sgn_my);
      if (sgn_mz!=0)
	sgn_mz/=fabs(sgn_mz);
      
      vector m1(sgn_mx,0.0,0.0);
      vector m1b(sgn_mx,0.0,0.0);
      vector m1f(sgn_mx,0.0,0.0); 
      for (int i=-1;i<2;i++)
	{
	  m1.y()-=(A[a3(i,1,0)]-A[a3(i,-1,0)])/2.0;
	  m1.z()-=(A[a3(i,0,1)]-A[a3(i,0,-1)])/2.0;
	}
      scalar m1S=fabs(m1.x())+fabs(m1.y())+fabs(m1.z());
	  
      vector m2(0.0,sgn_my,0.0);
      vector m2b(0.0,sgn_my,0.0);
      vector m2f(0.0,sgn_my,0.0); 
      for (int i=-1;i<2;i++)
	{
	  m2.x()-=(A[a3(1,i,0)]-A[a3(-1,i,0)])/2.0;
	  m2.z()-=(A[a3(0,i,1)]-A[a3(0,i,-1)])/2.0;
	}
      scalar m2S=fabs(m2.x())+fabs(m2.y())+fabs(m2.z());

      vector m3(0.0,0.0,sgn_mz);
      vector m3b(0.0,0.0,sgn_mz);
      vector m3f(0.0,0.0,sgn_mz); 
      for (int i=-1;i<2;i++)
	{
	  m3.x()-=(A[a3(1,0,i)]-A[a3(-1,0,i)])/2.0;
	  m3.y()-=(A[a3(0,1,i)]-A[a3(0,-1,i)])/2.0;
	}
      scalar m3S=fabs(m3.x())+fabs(m3.y())+fabs(m3.z());

      if (m1S==0)
	m1S=Foam::GREAT;
      if (m2S==0)
	m2S=Foam::GREAT;
      if (m3S==0)
	m3S=Foam::GREAT;

      
      if (fabs(m1.x())/m1S>fabs(m2.y())/m2S and fabs(m1.x())/m1S>fabs(m3.z())/m3S)
	{
	  return -m1;
	}
      else if (fabs(m2.y())/m2S>fabs(m1.x())/m1S and fabs(m2.y())/m2S>fabs(m3.z())/m3S)
	{
	  return -m2;
	}
      else
	{
	  return -m3;
	}
    }
}


void Foam::uniformStencil::flip(const int axis)
{
  List<scalar> A;
  A.setSize(N_,0.0);  
  
  int sgnnx=1;
  int sgnny=1;
  int sgnnz=1;
  if (axis==0)
    sgnnx=-1;
  else if(axis==1)
    sgnny=-1;
  else if(axis==2)
    sgnnz=-1;
  else
    assert(axis);

  for (int k=-kMax_;k<kMax_+1;k++)
    {
      for (int j=-jMax_;j<jMax_+1;j++)
	{
	  for (int i=-iMax_;i<iMax_+1;i++)
	    {
	      A[a3(i*sgnnx,j*sgnny,k*sgnnz)]=A_[a3(i,j,k)];
	    }
	}
    }
  A_=A;
}


void Foam::uniformStencil::transpose(const bool xaxis,const bool yaxis,const bool zaxis)
{
  assert(not(xaxis && yaxis && zaxis));

  List<scalar> A;
  A.setSize(N_,0.0);  

  if (xaxis && yaxis)
    {
      for (int k=-kMax_;k<kMax_+1;k++)
	for (int j=-jMax_;j<jMax_+1;j++)
	  for (int i=-iMax_;i<iMax_+1;i++)
	    A[a3(i,j,k)]=A_[a3(j,i,k)];
    }
  else if (xaxis && zaxis)
    {
      for (int k=-kMax_;k<kMax_+1;k++)
	for (int j=-jMax_;j<jMax_+1;j++)
	  for (int i=-iMax_;i<iMax_+1;i++)
	    A[a3(i,j,k)]=A_[a3(k,j,i)];
    }
  else if (yaxis && zaxis)
    {
      for (int k=-kMax_;k<kMax_+1;k++)
	for (int j=-jMax_;j<jMax_+1;j++)
	  for (int i=-iMax_;i<iMax_+1;i++)
	    A[a3(i,j,k)]=A_[a3(i,k,j)];
    }

  A_=A; 
}


Foam::List<Foam::scalar> Foam::uniformStencil::getSmallStencil(const int iMax,const int jMax,const int kMax,std::string format)
{
  Foam::List<scalar> A;

  if (format=="F")
    {
      for(int k=-kMax;k<kMax+1;k++)
	for(int j=-jMax;j<jMax+1;j++)
	  for(int i=-iMax;i<iMax+1;i++)
	    A.push_back(A_[a3(i,j,k)]);
    }
  else if (format=="C")
    {
      for(int i=-iMax;i<iMax+1;i++)
	for(int j=-jMax;j<jMax+1;j++)
	  for(int k=-kMax;k<kMax+1;k++)
	    A.push_back(A_[a3(i,j,k)]);
    }
  else
    assert("Flattening format unknown.");
      
  return A;
}

// ************************************************************************* //
