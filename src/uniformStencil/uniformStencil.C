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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

/*namespace Foam
{
  defineTypeNameAndDebug(uniformStencil, 0);
  }*/


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::uniformStencil::uniformStencil(const Foam::fvMesh& mesh,
				     const ijkZone& ijkMesh,
				     int iMax,
				     int jMax,
				     int kMax):
  mesh_(mesh),
  ijkMesh_(ijkMesh),
  is2D_(false),
  iMax_(iMax),
  jMax_(jMax),
  kMax_(kMax)
{  

  Info<<"uniformStencil: constructing..."<<endl;

  const Vector<bool>& isEmpty=ijkMesh_.isEmpty();
  if (isEmpty.x())
    {
      iMax_=0;is2D_=true;
    }
  if (isEmpty.y())
    {
      jMax_=0;is2D_=true;
    }
  if (isEmpty.z())
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

  //if (N_==NS_*NS_)
  //  is2D_=true;
  //if (iMax_==0 or jMax_==0 or kMax_==0)
  // is2D_=true;

}

Foam::uniformStencil::uniformStencil(const Foam::fvMesh& mesh,
				     const ijkZone& ijkMesh,
				     int nMax):
  mesh_(mesh),
  ijkMesh_(ijkMesh),
  //  nMax_(nMax),
  // NS_(2*nMax_+1),
  is2D_(false),
  iMax_(nMax),
  jMax_(nMax),
  kMax_(nMax),
  NxS_(2*nMax+1),
  NyS_(2*nMax+1),
  NzS_(2*nMax+1)
{  
  const Vector<bool>& isEmpty=ijkMesh_.isEmpty();
  if (isEmpty.x())
    {
      NxS_=1;is2D_=true;
    }
  if (isEmpty.y())
    {
      NyS_=1;is2D_=true;
    }
  if (isEmpty.z())
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
  const labelList& globalIds=ijkMesh_.globalIds();
  const Vector<bool>& isEmpty=ijkMesh_.isEmpty();
  const label Nx=ijkMesh_.Nx();
  const label Ny=ijkMesh_.Ny();
  const label Nz=ijkMesh_.Nz();
  //const bool isEmpty=ijkMesh_.symXIn();
  //const point cc = C[celli];
  
  const label& iP=ijk.x();
  const label& jP=ijk.y();
  const label& kP=ijk.z();

  int il,jl,kl=0;
      
  if (isEmpty.x())
    {
      il=0;
      for (int j=-jMax_;j<jMax_+1;j++)
	{
	  if ((jP+j+1>Ny and ijkMesh_.symYOut()) or (jP-j<0 and ijkMesh_.symYIn())) 
	    jl=jP-j;
	  else
	    jl=(jP+j+Ny)%Ny;
	  for (int k=-kMax_;k<kMax_+1;k++)
	    {
	      if ((kP+k+1>Nz and ijkMesh_.symZOut()) or (kP-k<0 and ijkMesh_.symZIn())) 
		kl=kP-k;
	      else
		kl=(kP+k+Nz)%Nz;
	      label ijkG=ijkMesh_.ijk1(il,jl,kl); 
	      label gblIdx=globalIds[ijkG];
	      A_[a2(j,k)]=phiIJK[gblIdx];	
	    }
	}
    }
  else if (isEmpty.y())
    {
      jl=0;
      for (int i=-iMax_;i<iMax_+1;i++)
	{
	  if ((iP+i+1>Nx and ijkMesh_.symXOut()) or (iP+i<0 and ijkMesh_.symXIn())) 
	    il=iP-i;
	  else
	    il=(iP+i+Nx)%Nx;
	  for (int k=-kMax_;k<kMax_+1;k++)
	    {
	      if ((kP+k+1>Nz and ijkMesh_.symZOut()) or (kP+k<0 and ijkMesh_.symZIn())) 
		kl=kP-k;
	      else
		kl=(kP+k+Nz)%Nz;
	      label ijkG=ijkMesh_.ijk1(il,jl,kl); 
	      label gblIdx=globalIds[ijkG];
	      A_[a2(i,k)]=phiIJK[gblIdx];
	    }
	}
    }
  else if (isEmpty.z())
    {
      kl=0;
      for (int i=-iMax_;i<iMax_+1;i++)
	{
	  if ((iP+i+1>Nx and ijkMesh_.symXOut()) or (iP+i<0 and ijkMesh_.symXIn())) 
	    il=iP-i;
	  else
	    il=(iP+i+Nx)%Nx;
	  for (int j=-jMax_;j<jMax_+1;j++)
	    {
	      if ((jP+j+1>Ny and ijkMesh_.symYOut()) or (jP+j<0 and ijkMesh_.symYIn())) 
		jl=jP-j;
	      else
		jl=(jP+j+Ny)%Ny;
	      label ijkG=ijkMesh_.ijk1(il,jl,kl); 
	      label gblIdx=globalIds[ijkG];
	      A_[a2(i,j)]=phiIJK[gblIdx];	
	    }
	}
    }
  else
    {
      for (int i=-iMax_;i<iMax_+1;i++)
	{
	  if ((iP+i+1>Nx and ijkMesh_.symXOut()) or (iP+i<0 and ijkMesh_.symXIn())) 
	    il=iP-i;
	  else
	    il=(iP+i+Nx)%Nx;
	  for (int j=-jMax_;j<jMax_+1;j++)
	    {
	      if ((jP+j+1>Ny and ijkMesh_.symYOut()) or (jP+j<0 and ijkMesh_.symYIn())) 
		jl=jP-j;
	      else
		jl=(jP+j+Ny)%Ny;
	      for (int k=-kMax_;k<kMax_+1;k++)
		{
		  if ((kP+k+1>Nz and ijkMesh_.symZOut()) or (kP+k<0 and ijkMesh_.symZIn())) 
		    kl=kP-k;
		  else
		    kl=(kP+k+Nz)%Nz;
		  //label il=(iP+i+Nx)%Nx;
		  //label jl=(jP+j+Ny)%Ny;
		  //label kl=(kP+k+Nz)%Nz;
		  label ijkG=ijkMesh_.ijk1(il,jl,kl); 
		  label gblIdx=globalIds[ijkG];
		  //label ijkL=(i+nMax_)+NS_*(j+nMax_)+NS_*NS_*(k+nMax_);
		  A_[a3(i,j,k)]=phiIJK[gblIdx];	
		}
	    }
	}
    }
}


Foam::scalar Foam::uniformStencil::estimateSignK()
{
  Foam::vector n=calcYoungNormal();
  return estimateSignK();
}


Foam::scalar Foam::uniformStencil::estimateSignK(Foam::vector n)
{
  const List<scalar>& A=A_;
  const Vector<bool>& isEmpty=ijkMesh_.isEmpty();
  scalar nN=0;
  scalar nT=0;
  label tMax=0;
  label nMax=0;
  scalar Hxx=0.0;
  if (is2D_)
    {      
      List<scalar> H;
      H.setSize(3,0.0);
      H=0;

      if (isEmpty.x())
	{
	  nT=n.y();tMax=jMax_;
	  nN=n.z();nMax=kMax_;
	}
      else if (isEmpty.y())
	{
	  nT=n.x();tMax=iMax_;
	  nN=n.z();nMax=kMax_;
	}
      else 
	{
	  nT=n.x();tMax=iMax_;
	  nN=n.y();nMax=jMax_;
	}
      
      if (Foam::mag(nN)>=Foam::mag(nT))
	{
	  for (int t=-1;t<2;t++)
	    {
	      for (int n=-nMax;n<nMax+1;n++)
		{
		  H[t+1]+=A[a2(t,n)];
		}
	    }
	}
      else
	{
	  for (int n=-1;n<2;n++)
	    {
	      for (int t=-tMax;t<tMax+1;t++)
		{
		  H[n+1]+=A[a2(t,n)];
		}
	    }
	}
      
      Hxx=(H[0]-2.*H[1]+H[2]);
    }
  else
    NotImplemented;

  return std::copysign(1.0,-Hxx); //(posK?1.0:-1.0);

}


Foam::vector Foam::uniformStencil::calcYoungNormal()
{
  const Foam::scalar dx=ijkMesh_.dx();
  Foam::vector m(0.0,0.0,0.0);

  const List<scalar>& A=A_;
  
  if (is2D_)
    {
      scalar mN=1./8./dx*(A[a2(1,1)]-A[a2(-1,1)]+2.0*A[a2(1,0)]-2.0*A[a2(-1,0)]+A[a2(1,-1)]-A[a2(-1,-1)]);
      scalar mT=1./8./dx*(A[a2(1,1)]-A[a2(1,-1)]+2.0*A[a2(0,1)]-2.0*A[a2(0,-1)]+A[a2(-1,1)]-A[a2(-1,-1)]);
      if (ijkMesh_.isEmpty().x())
	m=vector(0,mN,mT);
      else if (ijkMesh_.isEmpty().y())
	m=vector(mN,0,mT);
      else
	m=vector(mN,mT,0);
    }
  else
    {
      m.x()=1./32./dx*(A[a3(1,-1,-1)]-A[a3(-1,-1,-1)]+2*A[a3(1,0,-1)]-2*A[a3(-1,0,-1)]+A[a3(1,1,-1)]-A[a3(-1,1,-1)]+
		      2*A[a3(1,-1,0)]-2*A[a3(-1,-1,0)]+4*A[a3(1,0,0)]-4*A[a3(-1,0,0)]+2*A[a3(1,1,0)]-2*A[a3(-1,1,0)]+
		      A[a3(1,-1,1)]-A[a3(-1,-1,1)]+2*A[a3(1,0,1)]-2*A[a3(-1,0,1)]+A[a3(1,1,1)]-A[a3(-1,1,1)]);

      m.y()=1./32./dx*(A[a3(-1,1,-1)]-A[a3(-1,-1,-1)]+2*A[a3(0,1,-1)]-2*A[a3(0,-1,-1)]+A[a3(1,1,-1)]-A[a3(1,-1,-1)]+
		      2*A[a3(-1,1,0)]-2*A[a3(-1,-1,0)]+4*A[a3(0,1,0)]-4*A[a3(0,-1,0)]+2*A[a3(1,1,0)]-2*A[a3(1,-1,0)]+
		      A[a3(-1,1,1)]-A[a3(-1,-1,1)]+2*A[a3(0,1,1)]-2*A[a3(0,-1,1)]+A[a3(1,1,1)]-A[a3(1,-1,1)]);

      m.z()=1./32./dx*(A[a3(-1,-1,1)]-A[a3(-1,-1,-1)]+2*A[a3(-1,0,1)]-2*A[a3(-1,0,-1)]+A[a3(-1,1,1)]-A[a3(-1,1,-1)]+
		      2*A[a3(0,-1,1)]-2*A[a3(0,-1,-1)]+4*A[a3(0,0,1)]-4*A[a3(0,0,-1)]+2*A[a3(0,1,1)]-2*A[a3(0,1,-1)]+
		      A[a3(1,-1,1)]-A[a3(1,-1,-1)]+2*A[a3(1,0,1)]-2*A[a3(1,0,-1)]+A[a3(1,1,1)]-A[a3(1,1,-1)]);
    }

  //Info<<"calcYoungNormal: A= "<<A<<endl;
  //Info<<"calcYoungNormal: m= "<<m<<endl;
  
  return m;
}
 
Foam::vector Foam::uniformStencil::calcCCDNormal() 
{

  const List<scalar>& A=A_;

  if (is2D_)
    {
      Foam::vector m(0.0,0.0,0.0);
      scalar mxc,myc,mxb,myb,mxf,myf=0;
      scalar mx,my=0;
      scalar sgn_mx=0;
      scalar sgn_my=0;

      for (int i=-1;i<2;i++)
	{
	  sgn_mx-=A[a2(1,i)]-A[a2(-1,i)];
	  sgn_my-=A[a2(i,1)]-A[a2(i,-1)];
	}

      /*if(fabs(sgn_mx)==0)
	{
	  for (int i=-1;i<2;i++)
	    sgn_mx-=A[a2(0,i)]-A[a2(-1,i)];
	}
      if(fabs(sgn_my)==0)
	{
	  for (int i=-1;i<2;i++)
	    sgn_my-=A[a2(i,0)]-A[a2(i,-1)];
	    }*/

      if (sgn_mx!=0)
	sgn_mx/=fabs(sgn_mx);
      if (sgn_my!=0)
	sgn_my/=fabs(sgn_my);

      //sgn_mx/=fabs(sgn_mx);
      // sgn_my/=fabs(sgn_my);

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
      //if (fabs(mxf)>fabs(mx))
      //mx=mxf;
      //if (fabs(mxb)>fabs(mx))
      //mx=mxb;

      my=myc;
      //if (fabs(myf)>fabs(my))
      //my=myf;
      //if (fabs(myb)>fabs(my))
      //	my=myb;
      
      if (fabs(my)<fabs(mx))
	{
	  mx=sgn_mx;
	}
      else
	{
	  my=sgn_my;
	}


      if (ijkMesh_.isEmpty().x())
	m=vector(0,-mx,-my);
      else if (ijkMesh_.isEmpty().y())
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
	      //sgn_mxb-=A[a3(0,i,j)]-A[a3(-1,i,j)];
	      //sgn_myb-=A[a3(i,0,j)]-A[a3(i,-1,j)];
	      //sgn_mzb-=A[a3(i,j,0)]-A[a3(i,j,-1)];
	    }
	}

      //Info<<"signs: "<<sgn_mx<<sgn_my<<sgn_mz<<endl;

      /*if (sgn_mx==0)
	sgn_mx=sgn_mxb;
      if (sgn_my!=0)
	sgn_my=sgn_myb;
      if (sgn_mz==0)
      sgn_mz=sgn_mzb;*/

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
	  /*m1b.y()-=(A[a3(i,0,0)]-A[a3(i,-1,0)]);
	  m1b.z()-=(A[a3(i,0,0)]-A[a3(i,0,-1)]);
	  m1f.y()-=(A[a3(i,1,0)]-A[a3(i,0,0)]);
	  m1f.z()-=(A[a3(i,0,1)]-A[a3(i,0,0)]);*/
	}
      //if (fabs(m1b.y())+fabs(m1b.z())>fabs(m1.y())+fabs(m1.z()))
      //	  m1=m1b;
      //if (fabs(m1f.y())+fabs(m1f.z())>fabs(m1.y())+fabs(m1.z()))
      //	  m1=m1f;
      scalar m1S=fabs(m1.x())+fabs(m1.y())+fabs(m1.z());
	  
      vector m2(0.0,sgn_my,0.0);
      vector m2b(0.0,sgn_my,0.0);
      vector m2f(0.0,sgn_my,0.0); 
      for (int i=-1;i<2;i++)
	{
	  m2.x()-=(A[a3(1,i,0)]-A[a3(-1,i,0)])/2.0;
	  m2.z()-=(A[a3(0,i,1)]-A[a3(0,i,-1)])/2.0;
	  /*m2b.x()-=(A[a3(0,i,0)]-A[a3(-1,i,0)]);
	  m2b.z()-=(A[a3(0,i,0)]-A[a3(0,i,-1)]);
	  m2f.x()-=(A[a3(1,i,0)]-A[a3(0,i,0)]);
	  m2f.z()-=(A[a3(0,i,1)]-A[a3(0,i,0)]);*/
	}
      //if (fabs(m2b.x())+fabs(m2b.z())>fabs(m2.x())+fabs(m2.z()))
      //	  m2=m2b;
      //if (fabs(m2f.x())+fabs(m2f.z())>fabs(m2.x())+fabs(m2.z()))
      //	  m2=m2f;
      scalar m2S=fabs(m2.x())+fabs(m2.y())+fabs(m2.z());

      vector m3(0.0,0.0,sgn_mz);
      vector m3b(0.0,0.0,sgn_mz);
      vector m3f(0.0,0.0,sgn_mz); 
      for (int i=-1;i<2;i++)
	{
	  m3.x()-=(A[a3(1,0,i)]-A[a3(-1,0,i)])/2.0;
	  m3.y()-=(A[a3(0,1,i)]-A[a3(0,-1,i)])/2.0;
	  /*m3b.x()-=(A[a3(0,0,i)]-A[a3(-1,0,i)]);
	  m3b.y()-=(A[a3(0,0,i)]-A[a3(0,-1,i)]);
	  m3f.x()-=(A[a3(0,0,i)]-A[a3(0,0,i)]);
	  m3f.y()-=(A[a3(0,0,i)]-A[a3(0,0,i)]);*/
	}
      //if (fabs(m3b.x())+fabs(m3b.y())>fabs(m3.x())+fabs(m3.y()))
      //	  m3=m3b;
      //if (fabs(m3f.x())+fabs(m3f.y())>fabs(m3.x())+fabs(m3.y()))
      //	  m3=m3f;
      scalar m3S=fabs(m3.x())+fabs(m3.y())+fabs(m3.z());

      if (m1S==0)
	m1S=Foam::GREAT;
      if (m2S==0)
	m2S=Foam::GREAT;
      if (m3S==0)
	m3S=Foam::GREAT;

      
      if (fabs(m1.x())/m1S>fabs(m2.y())/m2S and fabs(m1.x())/m1S>fabs(m3.z())/m3S)
	{

	  //Info<<"Selected m1:  "<<"m1="<<m1/m1S<<",m2="<<m2/m2S<<",m3="<<m3/m3S<<",mY="<<mY/mYS<<endl;
	  //Info<<"Selected m1="<<-m1/m1S<<" ,mY="<<mY/mYS<<endl;
	  return -m1;
	}
      else if (fabs(m2.y())/m2S>fabs(m1.x())/m1S and fabs(m2.y())/m2S>fabs(m3.z())/m3S)
	{
	  //Info<<"Selected m2="<<-m2/m2S<<" ,mY="<<mY/mYS<<endl;
	  return -m2;
	}
      else
	{
	  //Info<<"Selected m3="<<-m3/m3S<<" ,mY="<<mY/mYS<<endl;
	  return -m3;
	}
    }
}


/*Foam::scalar Foam::uniformStencil::calcCurvature(label normalDir)
{
  const Foam::scalar dx=ijkMesh_.dx();
  scalar kappa=0.0;

  //const List<scalar>& A=A_;
  
  if (is2D_)
    {
      double H[3]={0.0,0.0,0.0}; 
      for (int iN=-nMax_;iN<nMax_+1;iN++)
	{
	  for (int iT=-1;iT<2;iT++)
	    {
	      if(normalDir==0)
		H[iT+1]+=A_[a2(iN,iT)]*dx;
	      else if(normalDir==1)
		{
		  if (ijkMesh_.isEmpty().x())
		    H[iT+1]+=A_[a2(iN,iT)]*dx;
		  else
		    H[iT+1]+=A_[a2(iT,iN)]*dx;
		}
	      else
		H[iT+1]+=A_[a2(iT,iN)]*dx;
	    }
	}
      scalar Ht=(H[2]-H[0])/2.0/dx;
      scalar Htt=(H[2]-2.0*H[1]+H[0])/dx/dx;
      kappa=(Htt)/Foam::pow(1.0+Ht*Ht,1.5);
    }
  else
    {
      double H[3][3]={{0.0,0.0,0.0},{0.0,0.0,0.0},{0.0,0.0,0.0}}; 
      for (int iN=-nMax_;iN<nMax_+1;iN++)//normal 
	{
	  for (int iT=-1;iT<2;iT++)//tangent
	    {
	      for (int iB=-1;iB<2;iB++)//bi-normal
		{
		  if (normalDir==0)
		    {
		      H[iT+1][iB+1]+=A_[a3(iN,iT,iB)]*dx;
		    }
		  else if (normalDir==1)
		    {
		      H[iT+1][iB+1]+=A_[a3(iT,iN,iB)]*dx;
		    }
		  else
		    {
		      H[iT+1][iB+1]+=A_[a3(iT,iB,iN)]*dx;
		    }			    			   
		}
	    }
	}
      scalar Ht=(H[2][1]-H[0][1])/2.0/dx;
      scalar Hb=(H[1][2]-H[1][0])/2.0/dx;
      scalar Htt=(H[2][1]-2.0*H[1][1]+H[0][1])/dx/dx;
      scalar Hbb=(H[1][2]-2.0*H[1][1]+H[1][0])/dx/dx;
      scalar Htb=(H[2][2]-H[2][0]-H[0][2]+H[0][0])/4./dx/dx;
      kappa=(Htt+Hbb+Htt*Hb*Hb+Hbb*Ht*Ht-2.0*Htb*Ht*Hb)/Foam::pow(1.0+Ht*Ht+Hb*Hb,1.5);      
    }

  return -kappa;
}
*/
/*std::vector<double> Foam::uniformStencil::getSmallStencil(int iMax, int jMax, int kMax)
{
  std::vector<double> data;
  for(int i=-iMax;i<iMax+1;i++)
    for(int j=-jMax;j<jMax+1;j++)
      for(int k=-kMax;k<kMax+1;k++)
	data.push_back(A_[a3(i,j,k)]);
  return data;
  }*/

// ************************************************************************* //
