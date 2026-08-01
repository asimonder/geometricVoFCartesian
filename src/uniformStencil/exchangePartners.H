// --- Exchange partners ---
        // Step 1: face-adjacent procs with sender cells (always correct)
        // Step 2: symmetric closure (prevents unmatched send/recv)
        // Step 3: tight bbox overlap for corner-adjacent procs
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

        // Symmetric closure via allNbrProcs
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

        // Tight bbox overlap for corner-adjacent procs
        // Uses potentialSenders bbox (tight in z) instead of domain bbox
        {
            point localSendMin = point::uniform(VGREAT);
            point localSendMax = point::uniform(-VGREAT);
            const pointField& cc = mesh_.C();
            forAll(potentialSenders, si)
            {
                const point& p = cc[potentialSenders[si]];
                localSendMin = Foam::min(localSendMin, p);
                localSendMax = Foam::max(localSendMax, p);
            }

            // Gather all procs' sender bboxes
            List<point> allSendMin(Pstream::nProcs(), point::uniform(VGREAT));
            List<point> allSendMax(Pstream::nProcs(), point::uniform(-VGREAT));
            allSendMin[Pstream::myProcNo()] = localSendMin;
            allSendMax[Pstream::myProcNo()] = localSendMax;
            Pstream::allGatherList(allSendMin);
            Pstream::allGatherList(allSendMax);

            const vector L(Lx_, Ly_, Lz_);

            for (label p = 0; p < Pstream::nProcs(); p++)
            {
                if (p == Pstream::myProcNo()) continue;
                if (exchangeSet.found(p)) continue;
                if (allSendMin[p].x() >= VGREAT/2) continue;

                // Check if their sender bbox overlaps with mine
                // (expanded by reach to cover stencil extent)
                point eMin = allSendMin[p] - point::uniform(reach);
                point eMax = allSendMax[p] + point::uniform(reach);

                if (bboxOverlaps(eMin, eMax,
                        localSendMin, localSendMax,
                        isCyclic_, L))
                {
                    exchangeSet.insert(p);
                }
            }
        }

        exchangeProcs_ = exchangeSet.toc();
        Foam::sort(exchangeProcs_);
