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

#include "uniformStencil.H"
//#include "dummyTransform.H"
//#include "emptyPolyPatch.H"
//#include "processorPolyPatch.H"
//#include "syncTools.H"
//#include "wedgePolyPatch.H"

//#include "globalPoints.H"


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
template<typename Type>
void Foam::uniformStencil::getZoneField(const Foam::boolList& zone,
				 Foam::Map<Type>& phiIJK,
				 const Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh>& phi,
				 Foam::Vector<label> stencilSize)
{

  List<labelHashSet> commNeeded(Pstream::nProcs());
  List<labelHashSet> send(Pstream::nProcs());
  const volVectorField& C = mesh_.C();
  label iMax=stencilSize.x();
  label jMax=stencilSize.y();
  label kMax=stencilSize.z();
  bool is2D=false;
  if (Ny_==1)
    {
      jMax=0;
      is2D=true;
    }
	      
  if (Nz_==1)
    {
      kMax=0;
      is2D=true;
    }

  forAll(zone,celli)
    {
      if (zone[celli])
	{
	  const vector cc = C[celli];
	  bool isActive=false;
	  if (Nz_==1)
	    {
	      if (cc.y()>(Pmin_.y()-1e-12+iMax*dy_) && cc.y()<(Pmax_.y()+1e-12-iMax*dy_))
		isActive=true;
	    }	      
	  else if (Ny_==1)
	    {
	      if (cc.z()>(Pmin_.z()-1e-12+iMax*dz_) && cc.z()<(Pmax_.z()+1e-12-iMax*dz_))
		isActive=true;
	    }
	  else
	    {
	      isActive=false;
	    }

	  label il,jl,kl=0;
	  if (isActive)
	    {
	      label i=round((cc.x()-Pmin_.x())/dx_);
	      label j=round((cc.y()-Pmin_.y())/dy_);
	      label k=round((cc.z()-Pmin_.z())/dz_);

	      for (int ii=-iMax;ii<iMax+1;ii++)
		{
		  for (int kk=-kMax;kk<kMax+1;kk++)
		    {
		      for (int jj=-jMax;jj<jMax+1;jj++)
			{
			  il=(i+ii+Nx_)%Nx_;
			  if (is2D)
			    {
			      jl=(j+jj);
			      kl=(k+kk); //make sure cellSet is big enough
			    }
			  else
			    {
			      jl=(j+jj+Ny_)%Ny_;
			      kl=(k+kk); //make sure cellSet is big enough
			    }
			  label ijk=il+Nx_*jl+Nx_*Ny_*kl;	  
			  label gblIdx=globalIds_[ijk];
			  
			  if (globalNumbering_.isLocal(gblIdx))
			    {
			      //Info<<"local Point: "<<cc<<endl;
			      if (!phiIJK.found(gblIdx))
				phiIJK.insert(gblIdx,
					      phi[globalNumbering_.toLocal(gblIdx)]);
			    }
			  else
			    {
			      //Info<<"non-local Point: "<<cc<<endl;
			      const label procID =
				globalNumbering_.whichProcID(gblIdx);
				  if (!commNeeded[procID].found(gblIdx))
				    commNeeded[procID].insert(gblIdx);
			    }
			}
		    }			  
		
		}
	    }
	}
    }

  //Info<<"getZoneField: done looping the interface... "<<endl;
  
  if (Pstream::parRun())
    {
      // Info<<"getZoneField: starting parallel communication... "<<endl;
      PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

      // Stream data into buffer
      for (label domain = 0; domain < Pstream::nProcs(); domain++)
	{
	  if (domain != Pstream::myProcNo())
	    {
	      // Put data into send buffer
	      UOPstream toDomain(domain, pBufs);
	      
	      toDomain << commNeeded[domain];
	    }
	}

  // wait until everything is written.
      pBufs.finishedSends();
  
      for (label domain = 0; domain < Pstream::nProcs(); domain++)
	{
	  send[domain].clear();
	  
	  if (domain != Pstream::myProcNo())
	    {
	      // get data from send buffer
	      UIPstream fromDomain(domain, pBufs);

	      fromDomain >> send[domain];
	    }
	}

      // Get values from other proc
      List<Map<Type>> sendValues(Pstream::nProcs());

      //Info<<"getZoneField: starting parallel communication: sendValues... "<<endl;
      forAll(send, domaini)
	{
	  for (const label sendIdx : send[domaini])
	    {
	      sendValues[domaini].insert
		(
		 sendIdx,
		 phi[globalNumbering_.toLocal(sendIdx)]
		 );
	    }
	}

      for (label domain = 0; domain < Pstream::nProcs(); domain++)
	{
	  if (domain != Pstream::myProcNo())
	    {
	      // Put data into send buffer
	      UOPstream toDomain(domain, pBufs);
	      
	      toDomain << sendValues[domain];
	    }
	}

      // wait until everything is written.
      pBufs.finishedSends();
      Map<Type> tmpValue;
      //Info<<"getZoneField: starting parallel communication: writing to phiIJK... "<<endl;
      for (label domain = 0; domain < Pstream::nProcs(); domain++)
	{
	  send[domain].clear();
	  
	  if (domain != Pstream::myProcNo())
	    {
	      // get data from send buffer
	      UIPstream fromDomain(domain, pBufs);

	      fromDomain >> tmpValue;
	      phiIJK+=tmpValue;
	    }
	}
  
    }

  //Info<<"getZoneField: done... "<<endl;
}


// ************************************************************************* //
