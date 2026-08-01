// commBboxDirect: face-adjacent + per-partner filtered exchange +
// round-2 relay for 2-hop diagonal coverage.
// Works with preservePatches.
//
// Partner detection: face-adjacent procs with symmetric closure.
//   O(nBoundaryFaces) detection, O(nFaceNeighbors) closure.
//   Catches procs with no interface that hold needed alpha data.
//
// Data exchange (two rounds):
//   Round 1 - per-partner cellInsideBBox filter over potentialSenders.
//     Each partner receives only cells within its interface reach.
//     Under preservePatches this filter is a no-op (interface bbox
//     spans the full domain in cyclic directions) but correct.
//   Round 2 - relay round-1 received data to face-adjacent partners
//     (excluding source). Covers 2-hop edge/diagonal procs.
//
// Potential senders: seeded from (a) interface zone, (b) cells
//   adjacent to processorPolyPatch, (c) cells adjacent to same-proc
//   cyclicPolyPatch. Expanded by expandLayers face-hops.
//
// Fills for commonTopologyFill.H:
//   spatialHash, remoteGIds, remoteCCs, gIdToValue, binSize, invBin,
//   phiMap, exchangeProcs_

{
    scalar tA0 = enableTimeLogs ? mesh_.time().elapsedCpuTime() : 0;

    // --- Estimate cell spacing from interfacial band ---
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

    const scalar reach = 1.5 * nLayers * dl;
    const vector L(Lx_, Ly_, Lz_);

    // --- Local interface bbox (for per-partner send filter) ---
    point localIMin = point::uniform(VGREAT);
    point localIMax = point::uniform(-VGREAT);
    forAll(zone, ci)
    {
        if (zone[ci])
        {
            const point& cc = mesh_.C()[ci];
            localIMin = Foam::min(localIMin, cc);
            localIMax = Foam::max(localIMax, cc);
        }
    }

    // --- Gather per-proc interface bboxes (for send filter) ---
    List<point> allIMin(Pstream::nProcs(), point::uniform(VGREAT));
    List<point> allIMax(Pstream::nProcs(), point::uniform(-VGREAT));
    allIMin[Pstream::myProcNo()] = localIMin;
    allIMax[Pstream::myProcNo()] = localIMax;
    Foam::reduce(allIMin, minOp<List<point>>());
    Foam::reduce(allIMax, maxOp<List<point>>());

    // --- Potential senders ---
    scalar tS0 = enableTimeLogs ? mesh_.time().elapsedCpuTime() : 0;

    bitSet isSender(mesh_.nCells(), false);
    forAll(zone, ci)
    {
        if (zone[ci]) isSender.set(ci);
    }

    {
        const polyBoundaryMesh& bMeshS = mesh_.boundaryMesh();
        forAll(bMeshS, patchi)
        {
            const polyPatch& pp = bMeshS[patchi];
            const bool isProc = isA<processorPolyPatch>(pp);
            const bool isCyc =
                isA<cyclicPolyPatch>(pp)
                && !isA<processorCyclicPolyPatch>(pp);
            if (!isProc && !isCyc) continue;

            const labelUList& fc = pp.faceCells();
            forAll(fc, fi)
            {
                isSender.set(fc[fi]);
            }
        }
    }

    const labelUList& own = mesh_.faceOwner();
    const labelUList& nei = mesh_.faceNeighbour();
    const label nInternal = mesh_.nInternalFaces();
    const label expandLayers =
        label(1.5 * max(max(iMax_, jMax_), kMax_)) + 2;

    for (label layer = 0; layer < expandLayers; layer++)
    {
        bitSet toAdd(mesh_.nCells(), false);
        for (label fi = 0; fi < nInternal; fi++)
        {
            const label o = own[fi];
            const label n = nei[fi];
            if (isSender.test(o) && !isSender.test(n))
                toAdd.set(n);
            else if (isSender.test(n) && !isSender.test(o))
                toAdd.set(o);
        }
        forAll(toAdd, ci)
        {
            if (toAdd.test(ci)) isSender.set(ci);
        }
    }

    labelList potentialSenders = isSender.toc();
    label nSenders = potentialSenders.size();
    const pointField& cc = mesh_.C();

    scalar tS1 = enableTimeLogs ? mesh_.time().elapsedCpuTime() : 0;
    scalar dtS = tS1 - tS0;

    // --- Face-adjacent partner detection + symmetric closure ---
    HashSet<label> exchangeSet;
    HashSet<label> allNbrProcs;
    {
        const polyBoundaryMesh& patches = mesh_.boundaryMesh();
        forAll(patches, patchi)
        {
            const polyPatch& pp = patches[patchi];
            if (!isA<processorPolyPatch>(pp)) continue;

            const processorPolyPatch& procPatch =
                refCast<const processorPolyPatch>(pp);
            const label nbrProc = procPatch.neighbProcNo();

            allNbrProcs.insert(nbrProc);

            if (exchangeSet.found(nbrProc)) continue;

            const labelUList& faceCells = pp.faceCells();
            forAll(faceCells, fi)
            {
                if (isSender.test(faceCells[fi]))
                {
                    exchangeSet.insert(nbrProc);
                    break;
                }
            }
        }
    }

    {
        PstreamBuffers pBufsSym(Pstream::commsTypes::nonBlocking);
        for (const label p : allNbrProcs)
        {
            UOPstream toNbr(p, pBufsSym);
            toNbr << exchangeSet.found(p);
        }
        pBufsSym.finishedSends();
        for (const label p : allNbrProcs)
        {
            bool nbrWantsMe;
            UIPstream fromNbr(p, pBufsSym);
            fromNbr >> nbrWantsMe;
            if (nbrWantsMe) exchangeSet.insert(p);
        }
    }

    exchangeProcs_ = exchangeSet.toc();
    Foam::sort(exchangeProcs_);

    scalar tA1 = enableTimeLogs ? mesh_.time().elapsedCpuTime() : 0;
    scalar dtA = tA1 - tA0;

    // --- Round 1: per-partner filtered exchange ---
    scalar tB0 = enableTimeLogs ? mesh_.time().elapsedCpuTime() : 0;

    DynamicList<label> allGIds;
    DynamicList<point> allCCs;

    labelList round1G;
    pointField round1C;
    List<Type> round1V;
    labelList round1Src;

    {
        PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

        forAll(exchangeProcs_, ni)
        {
            const label p = exchangeProcs_[ni];

            DynamicList<label> sendG;
            DynamicList<point> sendC;
            DynamicList<Type>  sendV;

            if (allIMin[p].x() < VGREAT/2)
            {
                const point eMin =
                    allIMin[p] - point::uniform(reach);
                const point eMax =
                    allIMax[p] + point::uniform(reach);

                forAll(potentialSenders, si)
                {
                    const label ci = potentialSenders[si];
                    if (cellInsideBBox(cc[ci], eMin, eMax,
                            isCyclic_, L))
                    {
                        sendG.append(globalNumbering_.toGlobal(ci));
                        sendC.append(cc[ci]);
                        sendV.append(phi[ci]);
                    }
                }
            }
            else
            {
                // Partner has no interface: send all potentialSenders
                // (partner may still hold cells needed via relay)
                forAll(potentialSenders, si)
                {
                    const label ci = potentialSenders[si];
                    sendG.append(globalNumbering_.toGlobal(ci));
                    sendC.append(cc[ci]);
                    sendV.append(phi[ci]);
                }
            }

            UOPstream toNbr(p, pBufs);
            toNbr << labelList(sendG)
                  << pointField(sendC)
                  << List<Type>(sendV);
        }

        pBufs.finishedSends();

        DynamicList<label> recvGAll;
        DynamicList<point> recvCAll;
        DynamicList<Type>  recvVAll;
        DynamicList<label> recvSrcAll;

        forAll(exchangeProcs_, ni)
        {
            labelList recvG;
            pointField recvC;
            List<Type> recvV;

            UIPstream fromNbr(exchangeProcs_[ni], pBufs);
            fromNbr >> recvG >> recvC >> recvV;

            forAll(recvG, ci)
            {
                const label gId = recvG[ci];
                if (!phiMap.found(gId))
                {
                    phiMap.insert(gId, recvV[ci]);
                }
                if (!gIdToValue.found(gId))
                {
                    gIdToValue.insert(gId, recvV[ci]);
                    allGIds.append(gId);
                    allCCs.append(recvC[ci]);
                    recvGAll.append(gId);
                    recvCAll.append(recvC[ci]);
                    recvVAll.append(recvV[ci]);
                    recvSrcAll.append(exchangeProcs_[ni]);
                }
            }
        }

        round1G = recvGAll;
        round1C = recvCAll;
        round1V = recvVAll;
        round1Src = recvSrcAll;
    }

    // --- Round 2: relay for 2-hop diagonal coverage ---
    {
        PstreamBuffers pBufs2(Pstream::commsTypes::nonBlocking);

        forAll(exchangeProcs_, ni)
        {
            const label p = exchangeProcs_[ni];

            DynamicList<label> sendG;
            DynamicList<point> sendC;
            DynamicList<Type>  sendV;

            forAll(round1G, ri)
            {
                if (round1Src[ri] != p)
                {
                    sendG.append(round1G[ri]);
                    sendC.append(round1C[ri]);
                    sendV.append(round1V[ri]);
                }
            }

            UOPstream toNbr(p, pBufs2);
            toNbr << labelList(sendG)
                  << pointField(sendC)
                  << List<Type>(sendV);
        }

        pBufs2.finishedSends();

        forAll(exchangeProcs_, ni)
        {
            labelList recvG;
            pointField recvC;
            List<Type> recvV;

            UIPstream fromNbr(exchangeProcs_[ni], pBufs2);
            fromNbr >> recvG >> recvC >> recvV;

            forAll(recvG, ci)
            {
                const label gId = recvG[ci];
                if (!phiMap.found(gId))
                {
                    phiMap.insert(gId, recvV[ci]);
                }
                if (!gIdToValue.found(gId))
                {
                    gIdToValue.insert(gId, recvV[ci]);
                    allGIds.append(gId);
                    allCCs.append(recvC[ci]);
                }
            }
        }
    }

    remoteGIds = allGIds;
    remoteCCs = allCCs;

    scalar tB1 = enableTimeLogs ? mesh_.time().elapsedCpuTime() : 0;
    scalar dtB = tB1 - tB0;

    // --- Spatial hash build ---
    scalar tC0 = enableTimeLogs ? mesh_.time().elapsedCpuTime() : 0;

    if (remoteCCs.size() > 0 && dl > SMALL)
    {
        binSize = dl;
        invBin = 1.0 / binSize;

        forAll(remoteCCs, ri)
        {
            const label key = packKey(spatialKey(remoteCCs[ri], invBin));
            if (!spatialHash.found(key))
            {
                spatialHash.insert(key, DynamicList<label>());
            }
            spatialHash[key].append(ri);
        }
    }

    scalar tC1 = enableTimeLogs ? mesh_.time().elapsedCpuTime() : 0;
    scalar dtC = tC1 - tC0;

    if (enableTimeLogs)
    {
        label nExch = exchangeProcs_.size();
        label nRemote = remoteCCs.size();

        reduce(dtA, maxOp<scalar>());
        reduce(dtS, maxOp<scalar>());
        reduce(dtB, maxOp<scalar>());
        reduce(dtC, maxOp<scalar>());
        reduce(nExch, maxOp<label>());
        reduce(nSenders, maxOp<label>());
        reduce(nRemote, maxOp<label>());

        Info<< "amrStencil[bboxDirect]: A(partners)=" << dtA
            << "s S(senders)=" << dtS
            << "s B(exchange)=" << dtB
            << "s C(spatHash)=" << dtC
            << "s | nExch=" << nExch
            << " nSenders=" << nSenders
            << " nRemote=" << nRemote << endl;
    }
}
