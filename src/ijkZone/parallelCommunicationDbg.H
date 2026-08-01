if (Pstream::parRun())
  {
    #include "calculateProcDistance.H"

    Map<label> procGlobalIds;
    if (cellIndices.size() > 0)
    {
      //globalIds_.setSize(Nx_ * Ny_ * Nz_, -1);       
      forAll(cellIndices, i)
       {
	 const label iCell = cellIndices[i];
	 label ijk = ijk1(ijk3(iCell));
	 label gId = globalNumbering_.toGlobal(iCell);
	 procGlobalIds.insert(ijk, gId);
       }
    }

    // Perform the two-way exchange with each neighbor.
    PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

    // Send your local map (which might be empty) to all your neighbors.
    forAll(neighbourProcs_, i)
      {
	UOPstream toNeighbor(neighbourProcs_[i], pBufs);
	toNeighbor << procGlobalIds;
      }
  
    // Wait for all non-blocking sends to complete before proceeding.
    pBufs.finishedSends();

    
    // Prepare to receive a map from each neighbor.
    List<Map<label>> receivedGlobalIDMaps(neighbourProcs_.size());
    forAll(neighbourProcs_, i)
      {
	UIPstream fromNeighbor(neighbourProcs_[i], pBufs);
	fromNeighbor >> receivedGlobalIDMaps[i];
      }

    forAllConstIter(Map<label>, procGlobalIds, iter)
      {
	PminL_ = point(VGREAT, VGREAT, VGREAT);
	PmaxL_ = point(-VGREAT, -VGREAT, -VGREAT);
	
	if(iter()==-1)
	  {
	    Pout<< "Proc " << Pstream::myProcNo()
		<< " -1 cell id detected!!!\n" << endl;
	  }
	const label lid = globalNumbering_.toLocal(iter());
	const point& p = C_[lid];
	PminL_ = Foam::min(PminL_, p);  
	PmaxL_ = Foam::max(PmaxL_, p);
      }
 
    if (cellIndices.size() > 0)
      {
	forAll(receivedGlobalIDMaps, i)
	  {
	    const Map<label>& receivedGidMap = receivedGlobalIDMaps[i];
	    forAllConstIter(Map<label>, receivedGidMap, iter)
	      {
		if(iter()==-1)
		  {
		  Pout<< "Proc " << Pstream::myProcNo()
		      << " -1 cell id detected!!!\n" << endl;
		  }
		label ijk=ijkGlobalToLocal(iter.key());
		Vector<label> ijk3=ijk1To3(ijk);
		const point p{dx_ * ijk3.x(), dy_ * ijk3.y(), dz_ * ijk3.z()};
		PminL_ = Foam::min(PminL_, p);  
		PmaxL_ = Foam::max(PmaxL_, p);
	      }
	  }
      }

    scalar Lx=PmaxL_.x()-PminL_.x();
    scalar Ly=PmaxL_.y()-PminL_.y();
    scalar Lz=PmaxL_.z()-PminL_.z();
    NxL_=round(Lx/dx_+1.0);
    NyL_=round(Ly/dy_+1.0);
    NzL_=round(Lz/dz_+1.0);

    if (cellIndices.size() > 0)
      {
	globalIds_.setSize(NxL_ * NyL_ * NzL_, -1);
      }
    
    forAllConstIter(Map<label>, procGlobalIds, iter)
      {
	if(iter()==-1)
	  {
	    Pout<< "Proc " << Pstream::myProcNo()
		<< " -1 cell id detected!!!\n" << endl;
	  }

	label ijk=ijkGlobalToLocal(iter.key());
	if (ijk < 0 || ijk >= NxL_*NyL_*NzL_) {
	  FatalErrorInFunction
	    << "Out of range local index: " << ijk
	    << " size: " << NxL_*NyL_*NzL_
	    << " proc: " << Pstream::myProcNo()
	    << exit(FatalError);
	}
	
	globalIds_[ijk] = iter();
      }
 
    if (cellIndices.size() > 0)
      {
	forAll(receivedGlobalIDMaps, i)
	  {
	    const Map<label>& receivedGidMap = receivedGlobalIDMaps[i];
	    forAllConstIter(Map<label>, receivedGidMap, iter)
	      {
		if(iter()==-1)
		  {
		  Pout<< "Proc " << Pstream::myProcNo()
		      << " -1 cell id detected!!!\n" << endl;
		  }
		label ijk=ijkGlobalToLocal(iter.key());
		if (ijk < 0 || ijk >= NxL_*NyL_*NzL_) {
		  FatalErrorInFunction
		    << "Out of range local index: " << ijk
		    << " size: " << NxL_*NyL_*NzL_
		    << " proc: " << Pstream::myProcNo()
		    << exit(FatalError);
		}

		globalIds_[ijk] = iter();
	      }
	  }
      }
    
  }

 else // This is the serial case
   {
    globalIds_.setSize(Nx_ * Ny_ * Nz_, -1);
    if (cellIndices.size() > 0)
    {
        forAll(cellIndices, i)
        {
            const label iCell = cellIndices[i];
            label ijk = ijk1(ijk3(iCell));
            globalIds_[ijk] = iCell;
        }
    }
}
