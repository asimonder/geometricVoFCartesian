  HashSet<label> neighbourProcs;
  if (cellIndices.size()>0)
    {
      const polyBoundaryMesh& boundary = mesh.boundaryMesh();

      forAll(boundary, patchi)
	{
	  const polyPatch& patch = boundary[patchi];
	
	  if (isA<processorPolyPatch>(patch))
	    {
	      const processorPolyPatch& procPatch = refCast<const processorPolyPatch>(patch);
	      label nbrProcNo = procPatch.neighbProcNo();
	      neighbourProcs.insert(nbrProcNo);
	    }
	}
    }

  if (Pstream::parRun())
    {
      if (cellIndices.size()>0)
 	{
	  //no 64 bit label to save memory
	  //List<int32_t> globalIds(Nx_*Ny_*Nz_,-1);
	  globalIds_.setSize(Nx_ * Ny_ * Nz_,-1);
      
	  //global ids on local processor
	  Map<label> procGlobalIds;
	  forAll(cellIndices, i)
	    {
	      const label iCell = cellIndices[i];
	      label ijk = ijk1(ijk3(iCell));
	      globalIds_[ijk]=globalNumbering_.toGlobal(iCell);
	      procGlobalIds.insert(ijk,globalIds_[ijk]);
	    }

	  // 4. Perform a two-way exchange with each neighbour.
	  PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);
	  
	  // --- SEND DATA ---
	  // Loop over your neighbours and send your local maps to each one.
	  forAll(neighbourProcs, i)
	    {
	      const label neighbourProcID = neighbourProcs[i];
	      UOPstream toNeighbour(neighbourProcID, pBufs);
	      // Send the map to the neighbour.
	      toNeighbour << procGlobalIds ;
	    }

	  // This command waits for all the non-blocking sends and receives to complete.
	  pBufs.finishedSends();
      
	  // --- RECEIVE DATA ---
	  // Create lists to hold the maps we will receive from our neighbours.
	  List<Map<label>> receivedGlobalIDMaps(neighbourProcs.size());
      
	  // Loop over your neighbours and post a non-blocking receive for each one.
	  forAll(neighbourProcs, i)
	    {
	      const label neighbourProcID = neighbourProcs[i];
	      UIPstream fromNeighbour(neighbourProcID, pBufs);
	      // Receive both maps from the neighbour.
	      fromNeighbour >> receivedGlobalIDMaps[i];
	    }
      

	  // 5. Merge the received maps into your main member variable maps.
	  //    This is a local operation with no communication.
	  forAll(receivedGlobalIDMaps, i)
	    {
	      const Map<label>& receivedGidMap = receivedGlobalIDMaps[i];
	      
	      forAllConstIter(Map<label>, receivedGidMap, it)
		{
		  const label ijk = it.key();
		  const label gId = it();
		  globalIds_[ijk] = gId;
		}
	    }
	}
    }
