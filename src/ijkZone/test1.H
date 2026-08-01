if (Pstream::parRun())
{
    // 1. Get the unique list of direct neighbors.
  const polyBoundaryMesh& boundary = mesh_.boundaryMesh();
  HashSet<label> neighbourSet;
  forAll(boundary, patchi)
    {
      const polyPatch& patch = boundary[patchi];
      if (isA<processorPolyPatch>(patch))
        {
          const processorPolyPatch& procPatch = refCast<const processorPolyPatch>(patch);
          neighbourSet.insert(procPatch.neighbProcNo());
        }
        else if (isA<processorCyclicPolyPatch>(patch))
        {
          const processorCyclicPolyPatch& cyclicPatch = refCast<const processorCyclicPolyPatch>(patch);
          neighbourSet.insert(cyclicPatch.neighbProcNo());
        }
    }
    const labelList neighbourProcs = neighbourSet.toc();


    ///////////// ALL PROCS///////////////////
    // const label myProc = Pstream::myProcNo();
    //const label nProcs = Pstream::nProcs();
    //labelList neighbourProcs(nProcs - 1);
    //label idx = 0;
    //for (label p = 0; p < nProcs; ++p)
    // {
    //if (p != myProc)
    //    {
    //      neighbourProcs[idx++] = p;
    //    }
    //}
    /////////////////////////////////////////

    Map<label> procGlobalIds;
    if (cellIndices.size() > 0)
    {
      // 4. Populate your final globalIds_ list.
      globalIds_.setSize(Nx_ * Ny_ * Nz_, -1); 

      forAll(cellIndices, i)
        {
            const label iCell = cellIndices[i];
            label ijk = ijk1(ijk3(iCell));
            label gId = globalNumbering_.toGlobal(iCell);
            procGlobalIds.insert(ijk, gId);
        }
    }

    // 3. Perform the two-way exchange with each neighbor.
    PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);
    // --- SEND DATA ---
    // Send your local map (which might be empty) to all your neighbors.
    forAll(neighbourProcs, i)
    {
        UOPstream toNeighbor(neighbourProcs[i], pBufs);
        toNeighbor << procGlobalIds;
    }

    // Wait for all non-blocking sends to complete before proceeding.
    pBufs.finishedSends();
    // --- RECEIVE DATA ---
    // Prepare to receive a map from each neighbor.

      List<Map<label>> receivedGlobalIDMaps(neighbourProcs.size());
      forAll(neighbourProcs, i)
	{
	  UIPstream fromNeighbor(neighbourProcs[i], pBufs);
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
