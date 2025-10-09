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

