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
      //      for (label domain = 0; domain < Pstream::nProcs(); domain++)
      forAll(neighbourProcs_, i)
	{
	  // Put data into send buffer
	  label domain=neighbourProcs_[i];
	  UOPstream toDomain(domain, pBufs);
	  toDomain << commNeeded[domain];
	}

  // wait until everything is written.
      pBufs.finishedSends();
  
      // for (label domain = 0; domain < Pstream::nProcs(); domain++)
      forAll(neighbourProcs_, i)
	{
	  label domain=neighbourProcs_[i];
	  send[domain].clear();
	  // get data from send buffer
	  UIPstream fromDomain(domain, pBufs);
	  fromDomain >> send[domain];
	}

      // Get values from other proc
      List<Map<Type>> sendValues(Pstream::nProcs());

      // Info<<"getZoneField: starting parallel communication: sendValues... "<<endl;
      //      forAll(send, domaini)
      forAll(neighbourProcs_, i)
	{
	  label domaini=neighbourProcs_[i];
	  for (const label sendIdx : send[domaini])
	    {
	      sendValues[domaini].insert
		(
		 sendIdx,
		 phi[globalNumbering_.toLocal(sendIdx)]
		 );
	    }
	}

      //      for (label domain = 0; domain < Pstream::nProcs(); domain++)
      forAll(neighbourProcs_, i)
	{
	      // Put data into send buffer
	  label	domain=neighbourProcs_[i];
	  UOPstream toDomain(domain, pBufs);
	  toDomain << sendValues[domain];
	}

      // wait until everything is written.
      pBufs.finishedSends();
      Map<Type> tmpValue;
      //Info<<"getZoneField: starting parallel communication: writing to phiIJK... "<<endl;
      //      for (label domain = 0; domain < Pstream::nProcs(); domain++)
      forAll(neighbourProcs_, i)
	{
	  label	domain=neighbourProcs_[i];
	  send[domain].clear();
	  
	  // get data from send buffer
	  UIPstream fromDomain(domain, pBufs);
	  
	  fromDomain >> tmpValue;
	  phiIJK+=tmpValue;
	}
  
    }

  //Info<<"getZoneField: done... "<<endl;
