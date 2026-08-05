// Face-lattice adjacency for the face-carried momentum predictor
// (DevPhiAdv). Built once for the static uniform lattice by hashing
// face centres; every internal axis-d face receives its same-axis
// neighbours along all axes and the four transverse advecting faces
// per transverse axis. Faces whose stencil is incomplete (physical
// boundaries, processor and cyclic seams, lattice edges) are marked
// notOk and keep the classical cell-interpolated predictor.

// lattice spacing per axis from the first internal face pair
vector hLat(vector::zero);
Vector<bool> axEmpty(false, false, false);
{
    scalar hMin =
        1.0/gMax(mesh.surfaceInterpolation::deltaCoeffs().primitiveField());
    for (direction d = 0; d < 3; ++d) { hLat[d] = hMin; }
}

// face axis, sign convention: faces grouped by dominant normal axis
const label nIntF = mesh.nInternalFaces();
labelList faceAxis(nIntF, -1);
{
    const surfaceVectorField& Sfl = mesh.Sf();
    const surfaceScalarField& magSfl = mesh.magSf();
    forAll(faceAxis, fi)
    {
        const vector nHat = Sfl[fi]/magSfl[fi];
        for (direction d = 0; d < 3; ++d)
        {
            if (mag(nHat[d]) > 0.99) { faceAxis[fi] = d; break; }
        }
    }
}

// an axis with no internal faces has a single cell layer: empty.
// This is exact, unlike any thickness heuristic.
{
    labelList nAx(3, label(0));
    forAll(faceAxis, fi)
    {
        if (faceAxis[fi] >= 0) { ++nAx[faceAxis[fi]]; }
    }
    for (direction d = 0; d < 3; ++d)
    {
        label n = nAx[d];
        reduce(n, sumOp<label>());
        axEmpty[d] = (n == 0);
    }
    Info<< "DevPhiAdv: empty axes = " << axEmpty << endl;
}

// hash of internal face centres: key = axis + packed lattice index
// of 2*Cf/h (face centres sit on half-integer lattice positions)
Map<label> faceHash;
{
    const surfaceVectorField& Cfl = mesh.Cf();
    const point bMin = boundBox(mesh.points(), true).min();
    for (label fi = 0; fi < nIntF; ++fi)
    {
        if (faceAxis[fi] < 0) continue;
        const point& c = Cfl[fi];
        const label ix = label(std::lround(2.0*(c.x()-bMin.x())/hLat.x()));
        const label iy = label(std::lround(2.0*(c.y()-bMin.y())/hLat.y()));
        const label iz = label(std::lround(2.0*(c.z()-bMin.z())/hLat.z()));
        const label key =
            faceAxis[fi] + 4*(ix + 100003*(iy + 100003*iz));
        faceHash.insert(key, fi);
    }
}

// adjacency arrays. For internal face fi of axis d:
//   fAdj[fi][2*e+0/1] : same-axis-d face shifted by -/+ h along axis e
//   uAdj[fi][4*e ...] : the four axis-e faces bracketing fi (for the
//                       transverse advecting velocity along e); for
//                       e == d the advecting velocity is phi itself
labelListList fAdj(nIntF);
labelListList uAdj(nIntF);
boolList faceOkAdv(nIntF, false);
{
    const surfaceVectorField& Cfl = mesh.Cf();
    const point bMin = boundBox(mesh.points(), true).min();

    auto lookup = [&](const point& c, const label axis) -> label
    {
        const label ix = label(std::lround(2.0*(c.x()-bMin.x())/hLat.x()));
        const label iy = label(std::lround(2.0*(c.y()-bMin.y())/hLat.y()));
        const label iz = label(std::lround(2.0*(c.z()-bMin.z())/hLat.z()));
        const label key = axis + 4*(ix + 100003*(iy + 100003*iz));
        const auto iter = faceHash.cfind(key);
        return iter.good() ? iter.val() : label(-1);
    };

    for (label fi = 0; fi < nIntF; ++fi)
    {
        const label d = faceAxis[fi];
        if (d < 0) continue;
        fAdj[fi].setSize(6, -1);
        uAdj[fi].setSize(12, -1);
        bool ok = true;
        const point& c = Cfl[fi];
        for (direction e = 0; e < 3; ++e)
        {
            if (axEmpty[e]) continue;      // empty axis contributes 0
            point pm = c; pm[e] -= hLat[e];
            point pp = c; pp[e] += hLat[e];
            fAdj[fi][2*e]   = lookup(pm, d);
            fAdj[fi][2*e+1] = lookup(pp, d);
            if (fAdj[fi][2*e] < 0 || fAdj[fi][2*e+1] < 0) { ok = false; }
            if (label(e) != d)
            {
                // four axis-e faces at c +/- h/2 in d and e
                label q = 0;
                for (label sd = -1; sd <= 1; sd += 2)
                {
                    for (label se = -1; se <= 1; se += 2)
                    {
                        point pq = c;
                        pq[d] += 0.5*sd*hLat[d];
                        pq[e] += 0.5*se*hLat[e];
                        uAdj[fi][4*e+q] = lookup(pq, e);
                        if (uAdj[fi][4*e+q] < 0) { ok = false; }
                        ++q;
                    }
                }
            }
        }
        faceOkAdv[fi] = ok;
    }

    label nOk = 0;
    forAll(faceOkAdv, fi) { if (faceOkAdv[fi]) ++nOk; }
    reduce(nOk, sumOp<label>());
    Info<< "DevPhiAdv: face-advection stencil complete on " << nOk
        << " internal faces; the rest keep the classical predictor"
        << endl;
}
