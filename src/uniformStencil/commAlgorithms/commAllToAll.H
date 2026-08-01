// allToAll: O(nProcs) bbox-based partner detection.
// Every proc checks bbox overlap with every other proc.
// Correct for all mesh types and cyclic BCs. Does NOT scale beyond
// a few hundred procs. Kept as trusted reference implementation.

{
    scalar tA0 = enableTimeLogs ? mesh_.time().elapsedCpuTime() : 0;

    // --- Gather interface bboxes ---
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

    List<point> allIMin(Pstream::nProcs(), point::uniform(VGREAT));
    List<point> allIMax(Pstream::nProcs(), point::uniform(-VGREAT));
    allIMin[Pstream::myProcNo()] = localIMin;
    allIMax[Pstream::myProcNo()] = localIMax;
    Foam::reduce(allIMin, minOp<List<point>>());
    Foam::reduce(allIMax, maxOp<List<point>>());

    // --- Gather domain bboxes ---
    const boundBox myBB = mesh_.bounds();
    List<point> allDMin(Pstream::nProcs(), point::uniform(VGREAT));
    List<point> allDMax(Pstream::nProcs(), point::uniform(-VGREAT));
    allDMin[Pstream::myProcNo()] = myBB.min();
    allDMax[Pstream::myProcNo()] = myBB.max();
    Foam::reduce(allDMin, minOp<List<point>>());
    Foam::reduce(allDMax, maxOp<List<point>>());

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

    const scalar reach = 1.5 * nLayers * dl;
    const vector L(Lx_, Ly_, Lz_);
    const bool iHaveInterface = (localIMin.x() < VGREAT/2);

    // --- Partner detection (symmetric) ---
    HashSet<label> exchangeSet;

    for (label p = 0; p < Pstream::nProcs(); p++)
    {
        if (p == Pstream::myProcNo()) continue;

        if (allIMin[p].x() < VGREAT/2)
        {
            point eMin = allIMin[p] - point::uniform(reach);
            point eMax = allIMax[p] + point::uniform(reach);
            if (bboxOverlaps(eMin, eMax, myBB.min(), myBB.max(),
                    isCyclic_, L))
            {
                exchangeSet.insert(p);
                continue;
            }
        }

        if (iHaveInterface)
        {
            point eMin = localIMin - point::uniform(reach);
            point eMax = localIMax + point::uniform(reach);
            if (bboxOverlaps(eMin, eMax, allDMin[p], allDMax[p],
                    isCyclic_, L))
            {
                exchangeSet.insert(p);
            }
        }
    }

    exchangeProcs_ = exchangeSet.toc();
    Foam::sort(exchangeProcs_);

    scalar tA1 = enableTimeLogs ? mesh_.time().elapsedCpuTime() : 0;
    scalar dtA = tA1 - tA0;

    // --- Per-partner exchange ---
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

            if (allIMin[p].x() < VGREAT/2)
            {
                point eMin = allIMin[p] - point::uniform(reach);
                point eMax = allIMax[p] + point::uniform(reach);

                for (label ci = 0; ci < mesh_.nCells(); ci++)
                {
                    if (cellInsideBBox(cc[ci], eMin, eMax,
                            isCyclic_, L))
                    {
                        sendG.append(globalNumbering_.toGlobal(ci));
                        sendC.append(cc[ci]);
                        sendV.append(phi[ci]);
                    }
                }
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
        label nRemote = remoteCCs.size();
        reduce(dtA, maxOp<scalar>());
        reduce(dtB, maxOp<scalar>());
        reduce(dtC, maxOp<scalar>());
        reduce(nExch, maxOp<label>());
        reduce(nRemote, maxOp<label>());

        Info<< "amrStencil[allToAll]: A(partners)=" << dtA
            << "s B(exchange)=" << dtB
            << "s C(spatHash)=" << dtC
            << "s | nExch=" << nExch
            << " nRemote=" << nRemote << endl;
    }
}
