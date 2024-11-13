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
//#include "dummyTransform.H"
//#include "emptyPolyPatch.H"
//#include "processorPolyPatch.H"
//#include "syncTools.H"
//#include "wedgePolyPatch.H"

//#include "globalPoints.H"


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
template<typename Type>
void Foam::ijkZone::getZoneField(const Foam::boolList& zone,
				 Foam::Map<Type>& phiIJK,
				 const Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh>& phi,
				 Foam::Vector<label> stencilSize)
{

  List<labelHashSet> commNeeded(Pstream::nProcs());
  List<labelHashSet> send(Pstream::nProcs());
  //const volVectorField& C = mesh_.C();
  label iMax=stencilSize.x();
  label jMax=stencilSize.y();
  label kMax=stencilSize.z();
  //bool is2D=false;

  if (isEmpty_.x())
    iMax=0;
  if (isEmpty_.y())
    jMax=0;
  if (isEmpty_.z())
    kMax=0;
  

  forAll(zone,celli)
    {
      if (zone[celli])
	{
	  //if (!atBoundary[celli])
	  Vector<label> ijkG=ijk3(celli);
	  label il,jl,kl=0;

	  label i=ijkG.x();
	  label j=ijkG.y();
	  label k=ijkG.z();

	  for (int ii=-iMax;ii<iMax+1;ii++)
	    {
	      for (int kk=-kMax;kk<kMax+1;kk++)
		{
		  for (int jj=-jMax;jj<jMax+1;jj++)
		    {
		      il=(i+ii+Nx_)%Nx_;
		      jl=(j+jj+Ny_)%Ny_;
		      kl=(k+kk+Nz_)%Nz_;
		      label ijk=il+Nx_*jl+Nx_*Ny_*kl;	  
		      label gblIdx=globalIds_[ijk];
		      
		      if (globalNumbering_.isLocal(gblIdx))
			{
			  //Info<<"local Point: "<<gblIdx<<endl;
			  //Info<<"phiIJK.size()="<<endl;
			  if (!phiIJK.found(gblIdx))
			    phiIJK.insert(gblIdx,
					  phi[globalNumbering_.toLocal(gblIdx)]);
			}
		      else
			{
			  //Info<<"non-local Point: "<<endl;
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

  //Info<<"getZoneField: done looping the interface... "<<endl;
  
  if (Pstream::parRun())
    {
      //Info<<"getZoneField: starting parallel communication... "<<endl;
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

      // Info<<"getZoneField: starting parallel communication: sendValues... "<<endl;
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
