// 1. Each processor calculates its LOCAL bounding box.
point localPmin = point::uniform(VGREAT);
point localPmax = point::uniform(-VGREAT);
if (cellIndices.size() > 0)
  {
    forAll(cellIndices, i)
      {
	const point& cc = C_[cellIndices[i]];
	localPmin = Foam::min(localPmin, cc);
	localPmax = Foam::max(localPmax, cc);
      }
  }
  
// Create lists to hold the data from all processors.
List<point> allPmins(Pstream::nProcs(), point::uniform(VGREAT)); // Default fill
List<point> allPmaxs(Pstream::nProcs(), point::uniform(-VGREAT));

// Put your local data into the lists.
allPmins[Pstream::myProcNo()] = localPmin;
allPmaxs[Pstream::myProcNo()] = localPmax;

Foam::reduce(allPmins, minOp<List<point>>());
Foam::reduce(allPmaxs, maxOp<List<point>>());

//Pstream::gatherList(allPmins);
//Pstream::broadcastList(allPmins);

// Then, gather the Pmax list.
//Pstream::gatherList(allPmaxs);
//Pstream::broadcastList(allPmaxs);

// Now, every processor has two lists: 'allPmins' and 'allPmaxs',
// containing the bounding box of every other processor.

// Define your communication threshold. This should be based on your
const scalar threshold = 10.0 * dx_; // Using the globally reduced dx_

// This will be your final list of processors to communicate with.
HashSet<label> procsToCommunicateWith;

if (cellIndices.size() > 0)
  {

    for (label procI = 0; procI < Pstream::nProcs(); ++procI)
      {
	// No need to check yourself or your direct neighbors again.
	if (procI == Pstream::myProcNo())
	  continue;
	
	if(allPmins[procI].x()>VGREAT/2. or allPmaxs[procI].x()<-VGREAT/2.)
	  continue;
	    
	// Get the bounding box of the remote processor.
	const point& remotePmin = allPmins[procI];
	const point& remotePmax = allPmaxs[procI];
	
	const vector L(Lx_, Ly_, Lz_);
	scalar sqrDist = 0.0;
    
	for (direction d = 0; d < vector::nComponents; ++d)
	  {
	    scalar delta = 0.0;
	    
	if (localPmax[d] < remotePmin[d]) // Case 1: A is to the left of B
	  {
	    // Direct distance (the gap)
	    delta = remotePmin[d] - localPmax[d];
	    
	    if (isCyclic_[d])
	      {
		// Correct wrapped distance calculation
		scalar span = remotePmax[d] - localPmin[d];
		scalar wrappedDelta = L[d] - span;
		delta = Foam::min(delta, wrappedDelta);
	      }
	  }
	else if (remotePmax[d] < localPmin[d]) // Case 2: B is to the left of A
	  {
	    // Direct distance (the gap)
	    delta = localPmin[d] - remotePmax[d];
	    
	    if (isCyclic_[d])
	      {
		// Correct wrapped distance calculation
		scalar span = localPmax[d] - remotePmin[d];
		scalar wrappedDelta = L[d] - span;
		delta = Foam::min(delta, wrappedDelta);
	      }
	  }
	// Case 3: Boxes overlap, delta remains 0.
	
	sqrDist += delta * delta;
	  }
	
	scalar dist = Foam::sqrt(sqrDist);      
	
	if (dist<threshold)
	  {
	    procsToCommunicateWith.insert(procI);
	  }
      }

  }
neighbourProcs_ = procsToCommunicateWith.toc();
