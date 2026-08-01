forAll(zone,celli)
{
    if (zone[celli])
    {
        Vector<label> ijkG = ijk3(celli);

        const label i = ijkG.x();
        const label j = ijkG.y();
        const label k = ijkG.z();

        bool skipCell = false;

        for (int ii = -iMax; ii <= iMax; ++ii)
        {
            for (int kk = -kMax; kk <= kMax; ++kk)
            {
                for (int jj = -jMax; jj <= jMax; ++jj)
                {
                    label iRaw = i + ii;
                    label jRaw = j + jj;
                    label kRaw = k + kk;

                    if (!isCyclic_.x() && !isEmpty_.x()
                        && !symXIn_ && iRaw < 0)
                        continue;
                    if (!isCyclic_.x() && !isEmpty_.x()
                        && !symXOut_ && iRaw >= Nx_)
                        continue;
                    if (!isCyclic_.y() && !isEmpty_.y()
                        && !symYIn_ && jRaw < 0)
                        continue;
                    if (!isCyclic_.y() && !isEmpty_.y()
                        && !symYOut_ && jRaw >= Ny_)
                        continue;
                    if (!isCyclic_.z() && !isEmpty_.z()
                        && !symZIn_ && kRaw < 0)
                        continue;
                    if (!isCyclic_.z() && !isEmpty_.z()
                        && !symZOut_ && kRaw >= Nz_)
                        continue;

                    label il = (iRaw + Nx_) % Nx_;
                    label jl = (jRaw + Ny_) % Ny_;
                    label kl = (kRaw + Nz_) % Nz_;

                    label ijk = ijkGlobalToLocal(il, jl, kl);

                    if (ijk < 0 || ijk >= globalIds_.size())
                        continue;

                    label gblIdx = globalIds_[ijk];

                    if (globalNumbering_.isLocal(gblIdx))
                    {
                        label lId = globalNumbering_.toLocal(gblIdx);

                        if (!phiIJK.found(gblIdx))
                        {
                            phiIJK.insert(gblIdx, phi[lId]);
                        }
                    }
                    else
                    {
                        const label procID = globalNumbering_.whichProcID(gblIdx);

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

