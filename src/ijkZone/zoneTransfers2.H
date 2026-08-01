
// This will hold the list of unique global IDs we need from other processors.
labelHashSet neededGlobalIds;

// This map will be populated with all the data we need for the stencil.
// We will pre-populate it with local values.
phiIJK.clear();

forAll(zone, celli)
{
    if (zone[celli])
    {
        Vector<label> ijkG = ijk3(celli);
        label i = ijkG.x();
        label j = ijkG.y();
        label k = ijkG.z();

        for (int ii = -iMax; ii < iMax + 1; ii++)
        {
            for (int kk = -kMax; kk < kMax + 1; kk++)
            {
                for (int jj = -jMax; jj < jMax + 1; jj++)
                {
                    label il = (i + ii + Nx_) % Nx_;
                    label jl = (j + jj + Ny_) % Ny_;
                    label kl = (k + kk + Nz_) % Nz_;
                    label ijk = il + Nx_*jl + Nx_*Ny_*kl;

                    // Look up the global ID from your pre-computed map.
                    label gblIdx = globalIds_[ijk];

                    // This check is critical to prevent the -1 error.
                    if (gblIdx == -1)
                    {
                        // This ijk is not part of the interface zone, so skip it.
                        continue;
                    }

                    // Check if the data is on the local processor.
                    if (globalNumbering_.isLocal(gblIdx))
                    {
                        // If it's local, get the value directly and add it to our final map.
                        if (!phiIJK.found(gblIdx))
                        {
                            phiIJK.insert(
                                gblIdx,
                                phi[globalNumbering_.toLocal(gblIdx)]
                            );
                        }
                    }
                    else
                    {
                        // If it's not local, add the global ID to the set of
                        // data we need to fetch. The HashSet automatically
                        // handles duplicates.
                        neededGlobalIds.insert(gblIdx);
                    }
                }
            }
        }
    }
}


// --- Step 2: Scalable Parallel Communication ---
if (Pstream::parRun())
  {

    // 2b. Sort the needed global IDs by which processor owns them.
    List<labelList> requestsToSend(Pstream::nProcs());
    const labelList neededGlobalIdsList = neededGlobalIds.toc();
    forAll(neededGlobalIdsList, i)
      {
	const label gId = neededGlobalIdsList[i];
	requestsToSend[globalNumbering_.whichProcID(gId)].append(gId);
      }

    // 2c. Perform the two-way, neighbor-only exchange.
    PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

    // --- PHASE 1: SEND REQUESTS ---
    // This part is correct. All processors send their request lists.
    forAll(neighbourProcs_, i)
      {
	const label neighborProcID = neighbourProcs_[i];
	UOPstream toNeighbor(neighborProcID, pBufs);
	toNeighbor << requestsToSend[neighborProcID];
      }

    pBufs.finishedSends();


    // --- PHASE 2: RECEIVE REQUESTS ---
    // ALL processors receive all pending requests from their neighbors.
    // We store them in a temporary list.
    List<labelList> receivedRequests(neighbourProcs_.size());
    forAll(neighbourProcs_, i)
      {
	const label neighborProcID = neighbourProcs_[i];
	UIPstream fromNeighbor(neighborProcID, pBufs);
	fromNeighbor >> receivedRequests[i];
      }
    
    
    // --- PHASE 3: SEND REPLIES ---
    // Now that all requests are received, we can fulfill them and send the data back.
    forAll(neighbourProcs_, i)
      {
	const label neighborProcID = neighbourProcs_[i];
	const labelList& requestedIds = receivedRequests[i];
	
	Map<Type> replyData;
	forAll(requestedIds, j)
	  {
	    const label gId = requestedIds[j];
	    replyData.insert(gId, phi[globalNumbering_.toLocal(gId)]);
	  }
	
	UOPstream toNeighbor(neighborProcID, pBufs);
	toNeighbor << replyData;
      }
    
    pBufs.finishedSends();


    // --- PHASE 4: RECEIVE REPLIES ---
    // Finally, receive the data you originally asked for.
    forAll(neighbourProcs_, i)
      {
	const label neighborProcID = neighbourProcs_[i];
	if (requestsToSend[neighborProcID].size() > 0)
	  {
	    UIPstream fromNeighbor(neighborProcID, pBufs);
	    Map<Type> receivedData;
	    fromNeighbor >> receivedData;
	    
	    // Merge the data into your final map.
	    for (typename Map<Type>::const_iterator iter = receivedData.begin();
		 iter != receivedData.end();
		 ++iter)
	      {
		phiIJK.insert(iter.key(), iter());
	      }
	  }
      }
  }










    
