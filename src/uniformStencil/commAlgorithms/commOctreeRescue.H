// commOctreeRescue: shared global octree over all cells.
// Each proc queries locally against the same tree, then fetches
// values from the owning proc in one round of exchange.

{
    scalar tE0 = enableTimeLogs ? mesh_.time().elapsedCpuTime() : 0;

    DynamicList<point> myMissingCC;
    DynamicList<label> myMissingCell;
    DynamicList<label> myMissingFlat;

    forAll(zone, celli)
    {
        if (!zone[celli]) continue;
        if (!stencilTopology_.found(celli)) continue;

        const labelList& topo = stencilTopology_[celli];
        const point seedCC = mesh_.C()[celli];
        const scalar V = mesh_.V()[celli];
        const scalar ldl = Foam::cbrt(V);

        for (label i = -iMax_; i <= iMax_; i++)
        for (label j = -jMax_; j <= jMax_; j++)
        for (label k = -kMax_; k <= kMax_; k++)
        {
            const label flatIdx = a3(i, j, k);
	    if (topo[flatIdx] >= 0) continue;

            point q = seedCC;
            if (!isEmpty_.x()) q.x() += i * ldl;
            if (!isEmpty_.y()) q.y() += j * ldl;
            if (!isEmpty_.z()) q.z() += k * ldl;

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

            myMissingCC.append(q);
            myMissingCell.append(celli);
            myMissingFlat.append(flatIdx);
        }
    }

    label nMissingGlobal = myMissingCC.size();
    reduce(nMissingGlobal, sumOp<label>());

    label nRescued = 0;

    if (nMissingGlobal > 0)
    {
        label localStale = (fineCellsVersion_ != mesh_.nCells()) ? 1 : 0;
        reduce(localStale, maxOp<label>());

        if (localStale == 1)
        {
            List<pointField> allCC(Pstream::nProcs());
            List<labelList>  allGId(Pstream::nProcs());

            pointField myCC(mesh_.nCells());
            labelList  myGId(mesh_.nCells());
            forAll(myCC, ci)
            {
                myCC[ci] = mesh_.C()[ci];
                myGId[ci] = globalNumbering_.toGlobal(ci);
            }
            allCC[Pstream::myProcNo()]  = myCC;
            allGId[Pstream::myProcNo()] = myGId;

            Pstream::allGatherList(allCC);
            Pstream::allGatherList(allGId);

            DynamicList<point> flatCC;
            DynamicList<label> flatGId;
            forAll(allCC, p)
            {
                forAll(allCC[p], i)
                {
                    flatCC.append(allCC[p][i]);
                    flatGId.append(allGId[p][i]);
                }
            }

            fineCellCC_      = flatCC;
            fineCellGId_     = flatGId;
            fineCellsVersion_ = mesh_.nCells();
        }

        scalar dl = Foam::cbrt(mesh_.V()[0]);
        reduce(dl, maxOp<scalar>());
        const scalar tolSq = Foam::sqr(0.3 * dl);

        Random rndGen(17301893);
        treeBoundBox bb(treeBoundBox(fineCellCC_).extend(rndGen, 1e-4));
        bb.min() -= point(1e-8, 1e-8, 1e-8);
        bb.max() += point(1e-8, 1e-8, 1e-8);

        indexedOctree<treeDataPoint> tree
        (
            treeDataPoint(fineCellCC_),
            bb, 8, 10, 3.0
        );

        List<DynamicList<label>> wantGId(Pstream::nProcs());
        List<DynamicList<label>> wantQi(Pstream::nProcs());

        forAll(myMissingCC, qi)
        {
            pointIndexHit hit = tree.findNearest(myMissingCC[qi], tolSq);
            if (!hit.hit()) continue;

            const label idx = hit.index();
            if (idx < 0) continue;

            const label gId = fineCellGId_[idx];
            const label owner = globalNumbering_.whichProcID(gId);

            if (owner == Pstream::myProcNo())
            {
                const label local = globalNumbering_.toLocal(gId);
                if (!phiMap.found(gId))
                {
                    phiMap.insert(gId, phi[local]);
                }
                stencilTopology_[myMissingCell[qi]][myMissingFlat[qi]] = gId;
                nRescued++;
            }
            else
            {
                wantGId[owner].append(gId);
                wantQi[owner].append(qi);
            }
        }

        PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

        for (label p = 0; p < Pstream::nProcs(); p++)
        {
            if (p == Pstream::myProcNo()) continue;

            UOPstream toNbr(p, pBufs);
            toNbr << labelList(wantGId[p]);
        }
        pBufs.finishedSends();

        PstreamBuffers pBufsReply(Pstream::commsTypes::nonBlocking);

        for (label p = 0; p < Pstream::nProcs(); p++)
        {
            if (p == Pstream::myProcNo()) continue;

            labelList reqGId;
            {
                UIPstream fromNbr(p, pBufs);
                fromNbr >> reqGId;
            }

            List<Type> reply(reqGId.size(), pTraits<Type>::zero);
            forAll(reqGId, i)
            {
                const label local = globalNumbering_.toLocal(reqGId[i]);
                reply[i] = phi[local];
            }

            UOPstream toNbr(p, pBufsReply);
            toNbr << reply;
        }
        pBufsReply.finishedSends();

        for (label p = 0; p < Pstream::nProcs(); p++)
        {
            if (p == Pstream::myProcNo()) continue;

            List<Type> reply;
            UIPstream fromNbr(p, pBufsReply);
            fromNbr >> reply;

            if (wantGId[p].empty()) continue;

            forAll(reply, i)
            {
                const label qi  = wantQi[p][i];
                const label gId = wantGId[p][i];

                if (!phiMap.found(gId))
                {
                    phiMap.insert(gId, reply[i]);
                }
                stencilTopology_[myMissingCell[qi]][myMissingFlat[qi]] = gId;
                nRescued++;
            }
        }
    }

    scalar tE1 = enableTimeLogs ? mesh_.time().elapsedCpuTime() : 0;
    scalar dtE = tE1 - tE0;
    if (enableTimeLogs)
    {
        reduce(dtE, maxOp<scalar>());
        reduce(nRescued, sumOp<label>());
        reduce(nBandCells, sumOp<label>());
        Info<< "amrStencil[rescue]: E(rescue)=" << dtE
            << "s | nBandCells=" << nBandCells
            << " nMissing=" << nMissingGlobal
            << " nRescued=" << nRescued << endl;
    }
}
