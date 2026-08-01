forAll(zone,celli)
{
    if (zone[celli])
    {
        Vector<label> ijkG = ijk3(celli);

        const label i = ijkG.x();
        const label j = ijkG.y();
        const label k = ijkG.z();

        for (int ii = -iMax; ii <= iMax; ++ii)
        {
            for (int kk = -kMax; kk <= kMax; ++kk)
            {
                for (int jj = -jMax; jj <= jMax; ++jj)
                {
                    label il = (i + ii + Nx_) % Nx_;
                    label jl = (j + jj + Ny_) % Ny_;
                    label kl = (k + kk + Nz_) % Nz_;

                    label ijk = ijkGlobalToLocal(il, jl, kl);

                    if (ijk < 0 || ijk >= NxL_ * NyL_ * NzL_)
                    {
                        FatalErrorInFunction
                            << "Out of range local index: " << ijk
                            << " size: " << NxL_ * NyL_ * NzL_
                            << " proc: " << Pstream::myProcNo()
                            << " (il,jl,kl)=(" << il << "," << jl << "," << kl << ")"
                            << exit(FatalError);
                    }

                    label gblIdx = globalIds_[ijk];

                    // Strict: 0 and -1 are invalid
                    if (gblIdx <= 0)
                    {
                        FatalErrorInFunction
                            << "Invalid global ID: " << gblIdx
                            << " at local index " << ijk
                            << " proc: " << Pstream::myProcNo()
                            << exit(FatalError);
                    }

                    if (globalNumbering_.isLocal(gblIdx))
                    {
                        label lId = globalNumbering_.toLocal(gblIdx);

                        if (lId < 0 || lId >= phi.size())
                        {
                            FatalErrorInFunction
                                << "Local id out of range: " << lId
                                << " phi.size()=" << phi.size()
                                << " for gblIdx=" << gblIdx
                                << " proc=" << Pstream::myProcNo()
                                << exit(FatalError);
                        }

                        if (!phiIJK.found(gblIdx))
                        {
                            phiIJK.insert(gblIdx, phi[lId]);
                        }
                    }
                    else
                    {
                        const label procID = globalNumbering_.whichProcID(gblIdx);

                        if (procID < 0 || procID >= commNeeded.size())
                        {
                            FatalErrorInFunction
                                << "Invalid procID=" << procID
                                << " for gblIdx=" << gblIdx
                                << " proc=" << Pstream::myProcNo()
                                << exit(FatalError);
                        }

                        if (!commNeeded[procID].found(gblIdx))
                        {
                            commNeeded[procID].insert(gblIdx);
                        }
                    }
                }
            }
        }
    }
}

/*forAll(zone,celli)
    {
      if (zone[celli])
	{
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
		      label ijk=ijkGlobalToLocal(il,jl,kl);//il+Nx_*jl+Nx_*Ny_*kl;
		      if ((ijk == 0) && (minIL_ > 0 || minJL_ > 0 || minKL_ > 0))
			{
			  FatalErrorInFunction
			    << "Proc " << Pstream::myProcNo()
			    << " got key=0 but local window starts at "
			    << "(" << minIL_ << "," << minJL_ << "," << minKL_ << ")"
			    << " => inconsistency in global/local mapping."
			    << exit(FatalError);
			}
		      if (ijk < 0 || ijk >= NxL_*NyL_*NzL_) {
			FatalErrorInFunction
			  << "Out of range local index: " << ijk
			  << " size: " << NxL_*NyL_*NzL_
			  << " proc: " << Pstream::myProcNo()
			  << exit(FatalError);
		      }

		      label gblIdx=globalIds_[ijk];
		      
		      if (globalNumbering_.isLocal(gblIdx))
			{
			  if (!phiIJK.found(gblIdx))
			    phiIJK.insert(gblIdx,
					  phi[globalNumbering_.toLocal(gblIdx)]);
			}
		      else
			{
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
*/


Info<<"getZoneField: done looping the interface... "<<endl;
  
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
