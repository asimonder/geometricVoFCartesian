// commRescueMissingStencils: exact-match rescue for stencil slots still
// at -1 after the main comm and commonTopologyFill. Orthogonal to the
// comm model (activated by rescueMissingStencils=true in amrStencilDict,
// regardless of which commAlgorithm is selected).
//
// Design:
// - Per-proc state only. No allGatherList of cell-level data, no global
//   octree. Memory O(nMissing + nIncomingQueries) per proc.
// - Exact-match tolerance. A query at position q matches a local cell ci
//   iff |cc[ci] - q| < exactTol (tight relative to dl). Near-matches are
//   rejected to avoid returning semantically wrong cells.
// - Two-round PstreamBuffers exchange: requests out, replies back.
// - Unmatched slots stay at -1, falling back to rdfCurvature via the
//   existing incomplete-stencil path.
//
// Consumes: stencilTopology_, phiMap, allIMin/allIMax from any comm
//   algorithm that gathered interface bboxes (all of them currently do,
//   via either local computation in the algo file or inherited state).
//   To stay independent of the algo, this file re-gathers the interface
//   bboxes cheaply (metadata only).

{
    scalar tR0 = enableTimeLogs ? mesh_.time().elapsedCpuTime() : 0;

    // --- Collect local missing queries ---
    // For each zone cell, for each stencil slot still at -1, compute the
    // cyclic-wrapped expected cell centre. Store (queryPoint, seedCelli,
    // flatIdx) so we can fill the slot once we get a reply.

    DynamicList<point> myQpos;
    DynamicList<label> myQseed;
    DynamicList<label> myQflat;

    forAll(zone, celli)
    {
        if (!zone[celli]) continue;
        if (!stencilTopology_.found(celli)) continue;

        const labelList& topo = stencilTopology_[celli];
        const point seedCC = mesh_.C()[celli];
        const scalar V = mesh_.V()[celli];

        scalar localDx, localDy, localDz;
        if (isEmpty_.x() && !isEmpty_.y() && !isEmpty_.z())
        {
            const scalar ld = Foam::sqrt(V / dx_);
            localDx = dx_; localDy = ld; localDz = ld;
        }
        else if (!isEmpty_.x() && isEmpty_.y() && !isEmpty_.z())
        {
            const scalar ld = Foam::sqrt(V / dy_);
            localDx = ld; localDy = dy_; localDz = ld;
        }
        else if (!isEmpty_.x() && !isEmpty_.y() && isEmpty_.z())
        {
            const scalar ld = Foam::sqrt(V / dz_);
            localDx = ld; localDy = ld; localDz = dz_;
        }
        else
        {
            const scalar ld = Foam::cbrt(V);
            localDx = ld; localDy = ld; localDz = ld;
        }

        for (label i = -iMax_; i <= iMax_; i++)
        for (label j = -jMax_; j <= jMax_; j++)
        for (label k = -kMax_; k <= kMax_; k++)
        {
            const label flatIdx = a3(i, j, k);
            if (topo[flatIdx] >= 0) continue;

            point q = seedCC;
            if (!isEmpty_.x()) q.x() += i * localDx;
            if (!isEmpty_.y()) q.y() += j * localDy;
            if (!isEmpty_.z()) q.z() += k * localDz;

            const point& bMin = domainMin_;
            if (isCyclic_.x() && Lx_ > SMALL)
            {
                if (q.x() < bMin.x()) q.x() += Lx_;
                else if (q.x() > bMin.x() + Lx_) q.x() -= Lx_;
            }
            if (isCyclic_.y() && Ly_ > SMALL)
            {
                if (q.y() < bMin.y()) q.y() += Ly_;
                else if (q.y() > bMin.y() + Ly_) q.y() -= Ly_;
            }
            if (isCyclic_.z() && Lz_ > SMALL)
            {
                if (q.z() < bMin.z()) q.z() += Lz_;
                else if (q.z() > bMin.z() + Lz_) q.z() -= Lz_;
            }

            myQpos.append(q);
            myQseed.append(celli);
            myQflat.append(flatIdx);
        }
    }

    label nMissingLocal = myQpos.size();
    label nMissingGlobal = nMissingLocal;
    reduce(nMissingGlobal, sumOp<label>());

    if (nMissingGlobal == 0)
    {
        if (enableTimeLogs)
        {
            scalar tR1 = mesh_.time().elapsedCpuTime();
            scalar dtR = tR1 - tR0;
            reduce(dtR, maxOp<scalar>());
            Info<< "amrStencil[rescueMissingStencils]: R=" << dtR
                << "s | nMissing=0 nRescued=0" << endl;
        }
    }
    else
    {
        // --- Estimate cell spacing for tolerance ---
        scalar dl = 0;
        forAll(zone, ci)
        {
            if (zone[ci])
            {
                const scalar V = mesh_.V()[ci];
                dl = (isEmpty_.x() || isEmpty_.y() || isEmpty_.z())
                    ? Foam::sqrt(V / max(dx_, max(dy_, dz_)))
                    : Foam::cbrt(V);
                break;
            }
        }
        reduce(dl, maxOp<scalar>());

        const scalar exactTol = 1e-6 * dl;
        const scalar exactTolSq = exactTol * exactTol;

        // --- Gather per-proc domain bboxes from local cell centres ---
        // Needed to decide which procs can possibly own a query position.
        point myDMin = point::uniform(VGREAT);
        point myDMax = point::uniform(-VGREAT);
        {
            const pointField& ccAll = mesh_.C();
            forAll(ccAll, ci)
            {
                myDMin = Foam::min(myDMin, ccAll[ci]);
                myDMax = Foam::max(myDMax, ccAll[ci]);
            }
        }

        List<point> allDMin(Pstream::nProcs(), point::uniform(VGREAT));
        List<point> allDMax(Pstream::nProcs(), point::uniform(-VGREAT));
        allDMin[Pstream::myProcNo()] = myDMin;
        allDMax[Pstream::myProcNo()] = myDMax;
        Foam::reduce(allDMin, minOp<List<point>>());
        Foam::reduce(allDMax, maxOp<List<point>>());

        // --- Build per-proc request lists ---
        // For each query, probe every proc's domain bbox (cyclic-aware).
        // Candidate procs are those whose bbox contains the query point
        // (or any cyclic image of it). Typically 1 candidate per query.
        const vector L(Lx_, Ly_, Lz_);

        List<DynamicList<point>> reqPos(Pstream::nProcs());
        List<DynamicList<label>> reqQid(Pstream::nProcs());

        forAll(myQpos, qi)
        {
            const point& q = myQpos[qi];
            for (label p = 0; p < Pstream::nProcs(); p++)
            {
                if (p == Pstream::myProcNo()) continue;
                if (cellInsideBBox(q, allDMin[p], allDMax[p],
                        isCyclic_, L))
                {
                    reqPos[p].append(q);
                    reqQid[p].append(qi);
                }
            }
        }

        // --- Round 1: send requests, receive requests from others ---
        List<pointField> incomingPos(Pstream::nProcs());
        {
            PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

            for (label p = 0; p < Pstream::nProcs(); p++)
            {
                if (p == Pstream::myProcNo()) continue;
                UOPstream toNbr(p, pBufs);
                toNbr << pointField(reqPos[p]);
            }
            pBufs.finishedSends();

            for (label p = 0; p < Pstream::nProcs(); p++)
            {
                if (p == Pstream::myProcNo()) continue;
                UIPstream fromNbr(p, pBufs);
                fromNbr >> incomingPos[p];
            }
        }

        // --- Local spatial hash over my own cells for exact-match lookup ---
        // Bin size = dl. Same bin-pack scheme as commonTopologyFill, so
        // the spatial keys match if the receiving proc had the cell.
        const scalar binLocal = dl;
        const scalar invBinLocal = 1.0 / binLocal;
        Map<DynamicList<label>> localHash;
        {
            const pointField& ccAll = mesh_.C();
            forAll(ccAll, ci)
            {
                const label key =
                    packKey(spatialKey(ccAll[ci], invBinLocal));
                if (!localHash.found(key))
                {
                    localHash.insert(key, DynamicList<label>());
                }
                localHash[key].append(ci);
            }
        }

        // --- Resolve incoming queries ---
        // For each query position from proc p, find the local cell whose
        // centre matches within exactTol. Reply with gId and value on
        // success, sentinel (gId = -1) on failure.
        List<labelList> replyGId(Pstream::nProcs());
        List<List<Type>> replyVal(Pstream::nProcs());

        const pointField& ccLocal = mesh_.C();

        for (label p = 0; p < Pstream::nProcs(); p++)
        {
            if (p == Pstream::myProcNo()) continue;

            const pointField& qs = incomingPos[p];
            labelList rG(qs.size(), -1);
            List<Type> rV(qs.size(), pTraits<Type>::zero);

            forAll(qs, qi)
            {
                const point& q = qs[qi];
                const Vector<label> ck = spatialKey(q, invBinLocal);

                label matchCell = -1;
                for (label di = -1; di <= 1 && matchCell < 0; di++)
                for (label dj = -1; dj <= 1 && matchCell < 0; dj++)
                for (label dk = -1; dk <= 1 && matchCell < 0; dk++)
                {
                    const label key = packKey
                    (
                        Vector<label>(ck.x()+di, ck.y()+dj, ck.z()+dk)
                    );
                    const auto iter = localHash.cfind(key);
                    if (!iter.good()) continue;

                    const DynamicList<label>& bucket = iter.val();
                    forAll(bucket, bi)
                    {
                        const label ci = bucket[bi];
                        if (Foam::magSqr(ccLocal[ci] - q) < exactTolSq)
                        {
                            matchCell = ci;
                            break;
                        }
                    }
                }

                if (matchCell >= 0)
                {
                    rG[qi] = globalNumbering_.toGlobal(matchCell);
                    rV[qi] = phi[matchCell];
                }
            }

            replyGId[p] = rG;
            replyVal[p] = rV;
        }

        // --- Round 2: send replies, receive replies ---
        List<labelList> inReplyGId(Pstream::nProcs());
        List<List<Type>> inReplyVal(Pstream::nProcs());
        {
            PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

            for (label p = 0; p < Pstream::nProcs(); p++)
            {
                if (p == Pstream::myProcNo()) continue;
                UOPstream toNbr(p, pBufs);
                toNbr << replyGId[p] << replyVal[p];
            }
            pBufs.finishedSends();

            for (label p = 0; p < Pstream::nProcs(); p++)
            {
                if (p == Pstream::myProcNo()) continue;
                UIPstream fromNbr(p, pBufs);
                fromNbr >> inReplyGId[p] >> inReplyVal[p];
            }
        }

        // --- Fill rescued slots ---
        label nRescuedLocal = 0;
        for (label p = 0; p < Pstream::nProcs(); p++)
        {
            if (p == Pstream::myProcNo()) continue;
            const labelList& rG = inReplyGId[p];
            const List<Type>& rV = inReplyVal[p];
            const DynamicList<label>& qids = reqQid[p];

            forAll(rG, i)
            {
                const label gId = rG[i];
                if (gId < 0) continue;

                const label qi = qids[i];
                const label celli = myQseed[qi];
                const label flatIdx = myQflat[qi];

                // Only fill if still unfilled: another proc may have
                // already resolved this slot (cyclic-image overlap).
                if (stencilTopology_[celli][flatIdx] >= 0) continue;

                stencilTopology_[celli][flatIdx] = gId;
                if (!phiMap.found(gId))
                {
                    phiMap.insert(gId, rV[i]);
                }
                nRescuedLocal++;
            }
        }

        if (enableTimeLogs)
        {
            scalar tR1 = mesh_.time().elapsedCpuTime();
            scalar dtR = tR1 - tR0;
            reduce(dtR, maxOp<scalar>());
            reduce(nMissingLocal, sumOp<label>());
            reduce(nRescuedLocal, sumOp<label>());
            Info<< "amrStencil[rescueMissingStencils]: R=" << dtR
                << "s | nMissing=" << nMissingLocal
                << " nRescued=" << nRescuedLocal << endl;
        }
    }
}
