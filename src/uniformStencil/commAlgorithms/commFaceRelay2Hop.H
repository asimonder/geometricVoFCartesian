// faceRelay2Hop: face-adjacent partners with symmetric closure,
// round 1 direct send, round 2 relay for 2-hop corner coverage.
// Covers face and edge adjacencies. Misses 3-hop corners (only
// relevant in full 3D decomposition with 8-proc meeting points).

{
    scalar tA0 = enableTimeLogs ? mesh_.time().elapsedCpuTime() : 0;

    bitSet isSender(mesh_.nCells(), false);
    forAll(zone, ci)
    {
        if (zone[ci]) isSender.set(ci);
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

    scalar tA1 = enableTimeLogs ? mesh_.time().elapsedCpuTime() : 0;
    scalar dtA = tA1 - tA0;

    scalar tB0 = enableTimeLogs ? mesh_.time().elapsedCpuTime() : 0;

    DynamicList<label> allGIds;
    DynamicList<point> allCCs;

    labelList round1G;
    pointField round1C;
    List<Type> round1V;
    labelList round1Src;

    // Round 1: send local potentialSenders
    {
        PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);
        const pointField& cc = mesh_.C();

        forAll(exchangeProcs_, ni)
        {
            const label p = exchangeProcs_[ni];

            DynamicList<label> sendG;
            DynamicList<point> sendC;
            DynamicList<Type>  sendV;

            forAll(potentialSenders, si)
            {
                const label ci = potentialSenders[si];
                sendG.append(globalNumbering_.toGlobal(ci));
                sendC.append(cc[ci]);
                sendV.append(phi[ci]);
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

    // Round 2: relay round-1 received data for 2-hop corner coverage
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
        label nSenders = potentialSenders.size();
        label nRemote = remoteCCs.size();
        reduce(dtA, maxOp<scalar>());
        reduce(dtB, maxOp<scalar>());
        reduce(dtC, maxOp<scalar>());
        reduce(nExch, maxOp<label>());
        reduce(nSenders, maxOp<label>());
        reduce(nRemote, maxOp<label>());

        Info<< "amrStencil[faceRelay2Hop]: A(senders+partners)=" << dtA
            << "s B(exchange)=" << dtB
            << "s C(spatHash)=" << dtC
            << "s | nExch=" << nExch
            << " nSenders=" << nSenders
            << " nRemote=" << nRemote << endl;
    }
}
