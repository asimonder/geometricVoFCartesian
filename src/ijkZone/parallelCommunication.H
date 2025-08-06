if (Pstream::parRun())
  {
    #include "calculateProcDistance.H"

    Map<label> procGlobalIds;
    if (cellIndices.size() > 0)
    {
      globalIds_.setSize(Nx_ * Ny_ * Nz_, -1);       
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
    
    label nLocal = 0;
    forAllConstIter(Map<label>, procGlobalIds, iter)
      {
	globalIds_[iter.key()] = iter();
	if(iter()==-1)
	  {
	    Pout<< "Proc " << Pstream::myProcNo()
		<< " -1 cell id detected!!!\n" << endl;
	  }
	++nLocal;
      }
 
    label nNonLocal = 0;
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
		globalIds_[iter.key()] = iter();
		++nNonLocal;
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
