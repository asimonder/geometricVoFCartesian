// faceOnly: face-adjacent partners with symmetric closure, single
// round of exchange. Misses corner neighbours reachable only via
// edge or point adjacency. Useful as baseline only.

{
    scalar tA0 = enableTimeLogs ? mesh_.time().elapsedCpuTime() : 0;

    // --- Build potentialSenders: zone expanded by expandLayers ---
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

    // --- Face-adjacent partner detection ---
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

    // --- Symmetric closure ---
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

    // --- Estimate cell spacing ---
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

    // --- Round 1 exchange ---
    scalar tB0 = enableTimeLogs ? mesh_.time().elapsedCpuTime() : 0;

    DynamicList<label> allGIds;
    DynamicList<point> allCCs;

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
                }
            }
        }
    }

    remoteGIds = allGIds;
    remoteCCs = allCCs;

    scalar tB1 = enableTimeLogs ? mesh_.time().elapsedCpuTime() : 0;
    scalar dtB = tB1 - tB0;

    // --- Build spatial hash ---
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

        Info<< "amrStencil[faceOnly]: A(senders+partners)=" << dtA
            << "s B(exchange)=" << dtB
            << "s C(spatHash)=" << dtC
            << "s | nExch=" << nExch
            << " nSenders=" << nSenders
            << " nRemote=" << nRemote << endl;
    }
}
