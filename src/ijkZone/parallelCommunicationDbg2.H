
if (Pstream::parRun())
{
    #include "calculateProcDistance.H"

    // key   = GLOBAL ijk (flattened with Nx_,Ny_)
    // value = GLOBAL cell id
    Map<label> procGlobalIds;

    // Collect local contributions
    if (cellIndices.size() > 0)
    {
        forAll(cellIndices, i)
        {
            const label iCell = cellIndices[i];                 // local id
            const label ijk  = ijk1(ijk3(iCell));               // GLOBAL ijk key
            const label gId  = globalNumbering_.toGlobal(iCell);// GLOBAL cell id
            procGlobalIds.insert(ijk, gId);
        }
    }

    // Exchange with neighbors
    PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

    forAll(neighbourProcs_, iSend)
    {
        UOPstream toNeighbor(neighbourProcs_[iSend], pBufs);
        toNeighbor << procGlobalIds;
    }
    pBufs.finishedSends();

    List<Map<label>> receivedGlobalIDMaps(neighbourProcs_.size());
    forAll(neighbourProcs_, iRecv)
    {
        UIPstream fromNeighbor(neighbourProcs_[iRecv], pBufs);
        fromNeighbor >> receivedGlobalIDMaps[iRecv];
    }

    // --------------------------
    // Build LOCAL integer window
    // --------------------------
    label minI = pTraits<label>::max, minJ = pTraits<label>::max, minK = pTraits<label>::max;
    label maxI = pTraits<label>::min, maxJ = pTraits<label>::min, maxK = pTraits<label>::min;

    // Local cells
    if (cellIndices.size() > 0)
    {
      forAll(cellIndices, ii)
	{
	  const label c = cellIndices[ii];
	  const Vector<label> g = ijk3(c);
	  minI = min(minI, g.x()); maxI = max(maxI, g.x());
	  minJ = min(minJ, g.y()); maxJ = max(maxJ, g.y());
	  minK = min(minK, g.z()); maxK = max(maxK, g.z());
	}
    }
      
    // Neighbors
    forAll(receivedGlobalIDMaps, r)
    {
        const Map<label>& m = receivedGlobalIDMaps[r];
        forAllConstIter(Map<label>, m, it)
        {
            const Vector<label> g = ijk1To3(it.key());
            minI = min(minI, g.x()); maxI = max(maxI, g.x());
            minJ = min(minJ, g.y()); maxJ = max(maxJ, g.y());
            minK = min(minK, g.z()); maxK = max(maxK, g.z());
        }
    }

    // Save integer offsets for local->global mapping
    minIL_ = minI;
    minJL_ = minJ;
    minKL_ = minK;

    maxIL_ = maxI;
    maxJL_ = maxJ;
    maxKL_ = maxK;

    
    const bool haveAnything = (minI != pTraits<label>::max);

    if (!haveAnything)
    {
      // no cells on this rank — just skip
        globalIds_.clear();
    }
    else
    {
        // Local window sizes (≥1)
        NxL_ = std::max<label>(1, maxI - minI + 1);
        NyL_ = std::max<label>(1, maxJ - minJ + 1);
        NzL_ = std::max<label>(1, maxK - minK + 1);

        // Align origins with global offset
        PminL_.x() = Pmin_.x() + dx_*minI;
        PminL_.y() = Pmin_.y() + dy_*minJ;
        PminL_.z() = Pmin_.z() + dz_*minK;

        PmaxL_.x() = Pmin_.x() + dx_*maxI;
        PmaxL_.y() = Pmin_.y() + dy_*maxJ;
        PmaxL_.z() = Pmin_.z() + dz_*maxK;

        // Allocate compact local table
        globalIds_.setSize(static_cast<size_t>(NxL_)*NyL_*NzL_, -1);

	
        // Fill with local entries
        forAllConstIter(Map<label>, procGlobalIds, it)
        {
            const label gId = it();          // GLOBAL cell id
            const label key = it.key();      // GLOBAL ijk key

	    const label ijkL = ijkGlobalToLocal(key);	    
	    
            globalIds_[ijkL] = gId;
        }

        // Fill with received entries
        forAll(receivedGlobalIDMaps, r2)
        {
            const Map<label>& m = receivedGlobalIDMaps[r2];
            forAllConstIter(Map<label>, m, it2)
            {
                const label gId = it2();
                const label key = it2.key();

                const label ijkL = ijkGlobalToLocal(key);
		
                globalIds_[ijkL] = gId;
            }
        }

	
    }
}
else // serial
{
    globalIds_.setSize(Nx_*Ny_*Nz_, -1);
    if (cellIndices.size() > 0)
    {
        forAll(cellIndices, i)
        {
            const label iCell = cellIndices[i];
            const label ijk  = ijk1(ijk3(iCell)); // GLOBAL==LOCAL in serial
            globalIds_[ijk]  = iCell;
        }
    }
}
