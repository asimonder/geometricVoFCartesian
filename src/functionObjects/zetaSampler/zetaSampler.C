/*---------------------------------------------------------------------------*\
            Copyright (c) 2026, Asim Onder
\*---------------------------------------------------------------------------*/

#include "zetaSampler.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "fvMesh.H"
#include "Pstream.H"
#include "PstreamReduceOps.H"
#include "OSspecific.H"
#include "interpolationCellPoint.H"
#include "Map.H"
#include "IOList.H"
#include "IOdictionary.H"
#include "boundBox.H"
#include <fstream>
#include <sstream>
#include <iomanip>
#include <limits>
#include <cmath>


namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(zetaSampler, 0);
    addToRunTimeSelectionTable(functionObject, zetaSampler, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::zetaSampler::zetaSampler
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    etaSourceName_("surfaceElevation"),
    dz0_(0),
    stretch_(0),
    zMaxAir_(0),
    zMaxWater_(0),
    kp_(0),
    xmin_(0),
    ymin_(0),
    Lx_(0),
    Ly_(0),
    flatInterface_(false),
    tiersAir_(),
    tiersWater_(),
    fieldNames_(),
    zetaLevels_(),
    levelNxy_(),
    levelDecay_(),
    outDirBase_(),
    zetaPlanesFileWritten_(false),
    fieldsChecked_(false)
{
    read(dict);

    outDirBase_ =
        mesh_.time().globalPath() / "postProcessing" / name;

    if (Pstream::master())
    {
        mkDir(outDirBase_);
    }

    Info<< "zetaSampler: configured" << nl
        << "  etaSource     = " << etaSourceName_ << nl
        << "  flatInterface = "
        << (flatInterface_
            ? "true (rigid-lid, z* = zeta_k)"
            : "false (wave-following, z* = zeta_k + eta*exp(-kp|zeta_k|))")
        << nl
        << "  kp            = " << kp_
        << (kp_ > 0
            ? " rad/m (Hara-Sullivan damped mapping, e-folding 1/kp)"
            : " (uniform shift, z* = eta + zeta_k)")
        << nl
        << "  xmin, ymin    = " << xmin_ << ", " << ymin_ << nl
        << "  Lx, Ly        = " << Lx_   << ", " << Ly_   << nl
        << "  dz0, stretch  = " << dz0_ << ", " << stretch_ << nl
        << "  zMaxAir/Water = " << zMaxAir_ << ", " << zMaxWater_ << nl
        << "  nTiersAir     = " << tiersAir_.size() << nl
        << "  nTiersWater   = " << tiersWater_.size() << nl
        << "  nZeta         = " << zetaLevels_.size() << nl
        << "  fields        = " << fieldNames_ << endl;

    {
        scalar prevZMax = 0;
        forAll(tiersAir_, t)
        {
            Info<< "  tierAir   " << t
                << " : Nxy=" << tiersAir_[t].Nxy
                << "  zeta in [" << prevZMax
                << ", " << tiersAir_[t].zMaxAbs << "]"
                << endl;
            prevZMax = tiersAir_[t].zMaxAbs;
        }
    }
    {
        scalar prevZMax = 0;
        forAll(tiersWater_, t)
        {
            Info<< "  tierWater " << t
                << " : Nxy=" << tiersWater_[t].Nxy
                << "  zeta in [" << -tiersWater_[t].zMaxAbs
                << ", " << -prevZMax << "]"
                << endl;
            prevZMax = tiersWater_[t].zMaxAbs;
        }
    }

    // Write the level distribution to disk so the user can plot it
    // immediately, without waiting for first execute().
    writeZetaLevelsFile();
}


// * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * * //

bool Foam::functionObjects::zetaSampler::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    dict.readIfPresent("etaSource",     etaSourceName_);
    dict.readEntry("dz0",               dz0_);
    dict.readEntry("stretch",           stretch_);
    dict.readEntry("zMaxAir",           zMaxAir_);
    dict.readEntry("zMaxWater",         zMaxWater_);
    // Decay rate of the wave-following deformation (Hara & Sullivan
    // 2015). The plane at level zeta_k deforms by D_k * eta with
    // D_k = exp(-kp * |zeta_k|), so planes follow the surface fully
    // at zeta = 0 and relax to flat fixed-z planes over the scale
    // 1/kp. Optional; kp = 0 (default) gives D_k = 1 everywhere and
    // recovers the uniform-shift mapping exactly. Monotonicity of the
    // mapping requires kp * max|eta| < 1; checked fatally on the
    // instantaneous eta grid at every execute.
    kp_ = 0;
    dict.readIfPresent("kp", kp_);
    // Optional: ignore eta(x, y, t) and treat the surface as a flat
    // plane at z = z0 (typically still water level). Sample targets
    // become z* = zeta_k instead of z* = zeta_k + D_k * eta. Useful
    // for a "rigid-lid" comparison where horizontal-plane statistics
    // make sense above the highest crest / below the lowest trough,
    // and for Favre-style horizontal averaging across the full domain.
    flatInterface_ = false;
    dict.readIfPresent("flatInterface", flatInterface_);
    // Horizontal extent comes from the mesh's global bounding box.
    // mesh_.bounds() returns the box already reduced across all
    // procs, so it gives the global domain extent.
    const boundBox& meshBB = mesh_.bounds();
    xmin_ = meshBB.min().x();
    ymin_ = meshBB.min().y();
    Lx_   = meshBB.max().x() - meshBB.min().x();
    Ly_   = meshBB.max().y() - meshBB.min().y();
    dict.readEntry("fields",            fieldNames_);

    if (dz0_ <= 0 || stretch_ < 1 || zMaxAir_ < 0 || zMaxWater_ < 0)
    {
        FatalIOErrorInFunction(dict)
            << "Invalid level distribution parameters: "
            << "need dz0 > 0, stretch >= 1, zMaxAir >= 0, zMaxWater >= 0."
            << exit(FatalIOError);
    }

    if (kp_ < 0)
    {
        FatalIOErrorInFunction(dict)
            << "Invalid decay rate: need kp >= 0 (kp = 0 recovers the "
            << "uniform shift)."
            << exit(FatalIOError);
    }

    // Parse two tier lists: tiersAir and tiersWater. Each entry is
    // (Nxy zMaxAbs). Both lists sorted ascending by zMaxAbs at read.
    auto parseTierList =
    [&](const Foam::word& key, Foam::List<Tier>& out)
    {
        out.clear();
        ITstream& is = dict.lookup(key);
        token tok(is);
        if (!tok.isPunctuation() || tok.pToken() != token::BEGIN_LIST)
        {
            FatalIOErrorInFunction(dict)
                << "Expected '(' to begin " << key << " list"
                << exit(FatalIOError);
        }

        DynamicList<Tier> tmp;
        while (true)
        {
            tok = token(is);
            if (tok.isPunctuation() && tok.pToken() == token::END_LIST) break;
            is.putBack(tok);

            Tier T;
            token open(is);
            if (!open.isPunctuation() || open.pToken() != token::BEGIN_LIST)
            {
                FatalIOErrorInFunction(dict)
                    << "Expected '(' to begin tier entry in " << key
                    << exit(FatalIOError);
            }
            is >> T.Nxy >> T.zMaxAbs;
            token close(is);
            if (!close.isPunctuation() || close.pToken() != token::END_LIST)
            {
                FatalIOErrorInFunction(dict)
                    << "Expected ')' to end tier entry in " << key
                    << exit(FatalIOError);
            }
            tmp.append(T);
        }

        Foam::sort
        (
            tmp,
            [](const Tier& a, const Tier& b){ return a.zMaxAbs < b.zMaxAbs; }
        );
        out = tmp;
    };
    parseTierList("tiersAir",   tiersAir_);
    parseTierList("tiersWater", tiersWater_);

    if (tiersAir_.empty() || tiersWater_.empty())
    {
        FatalIOErrorInFunction(dict)
            << "tiersAir and tiersWater must both be non-empty"
            << exit(FatalIOError);
    }

    buildZetaLevels();

    // Assign each level to a tier (uses signed zeta for air/water choice)
    levelNxy_.setSize(zetaLevels_.size(), -1);
    forAll(zetaLevels_, k)
    {
        const label N = tierNxy(zetaLevels_[k]);
        if (N <= 0)
        {
            FatalErrorInFunction
                << "zeta level " << zetaLevels_[k]
                << " falls outside every tier."
                << exit(FatalError);
        }
        levelNxy_[k] = N;
    }

    // Precompute the per-level decay factors D_k = exp(-kp * |zeta_k|)
    // of the Hara-Sullivan mapping. Constant per level (kp is a fixed
    // dictionary scalar), so they are evaluated once here and threaded
    // into the sampling loop and the manifest. kp = 0 gives all ones.
    levelDecay_.setSize(zetaLevels_.size(), 1.0);
    forAll(zetaLevels_, k)
    {
        levelDecay_[k] = Foam::exp(-kp_ * mag(zetaLevels_[k]));
    }

    // NOTE: field-existence check moved to first execute().
    // The grad functionObject (and similar derived-field producers)
    // does not register its result field in the mesh objectRegistry
    // until its own first execute(). At zetaSampler construction time
    // those fields therefore don't exist yet, even when the producing
    // functionObject is correctly listed earlier in controlDict.

    return true;
}


Foam::label Foam::functionObjects::zetaSampler::tierNxy
(
    scalar zeta
) const
{
    // Air side (zeta > 0) uses tiersAir_; water side uses tiersWater_.
    // For zeta = 0 we pick the air list (innermost tier on either side
    // is the same Nxy in well-formed configs).
    const List<Tier>& list = (zeta >= 0) ? tiersAir_ : tiersWater_;
    const scalar absZeta = mag(zeta);
    forAll(list, t)
    {
        if (absZeta <= list[t].zMaxAbs + SMALL)
            return list[t].Nxy;
    }
    return -1;
}


void Foam::functionObjects::zetaSampler::buildZetaLevels()
{
    DynamicList<scalar> tmp;

    // zeta = 0 level
    tmp.append(0.0);

    // Air side: walk +zeta with geometric spacing
    {
        scalar dz = dz0_;
        scalar z = 0;
        while (true)
        {
            z += dz;
            if (z > zMaxAir_ + SMALL) break;
            tmp.append(z);
            dz *= stretch_;
        }
    }

    // Water side: walk -zeta with geometric spacing
    {
        scalar dz = dz0_;
        scalar z = 0;
        while (true)
        {
            z += dz;
            if (z > zMaxWater_ + SMALL) break;
            tmp.append(-z);
            dz *= stretch_;
        }
    }

    Foam::sort(tmp);
    zetaLevels_ = scalarList(tmp);

    Info<< "zetaSampler: built " << zetaLevels_.size()
        << " zeta levels  (dz0=" << dz0_
        << " stretch=" << stretch_
        << " zMaxAir=" << zMaxAir_
        << " zMaxWater=" << zMaxWater_ << ")"
        << "  range=[" << zetaLevels_.first()
        << ", " << zetaLevels_.last() << "]" << endl;
}


void Foam::functionObjects::zetaSampler::writeZetaLevelsFile() const
{
    if (!Pstream::master()) return;

    std::ofstream os((outDirBase_ / "zetaLevels.txt").c_str());
    os << "# zeta levels, per-level Nxy and decay factor D_k\n"
       << "# Generated from dz0=" << dz0_
       << " stretch=" << stretch_
       << " zMaxAir=" << zMaxAir_
       << " zMaxWater=" << zMaxWater_
       << " kp=" << kp_ << "\n"
       << "# D_k = exp(-kp*|zeta_k|); plane k samples"
       << " z* = zeta_k + D_k*eta(x, y, t)\n"
       << "# columns: index   zeta [m]   Nxy   decay\n";
    forAll(zetaLevels_, k)
    {
        os << std::setw(6) << k
           << "   " << std::scientific << std::setprecision(8)
           << std::setw(16) << zetaLevels_[k]
           << "   " << std::setw(6) << levelNxy_[k]
           << "   " << std::setw(16) << levelDecay_[k]
           << "\n";
    }
    os.close();
}


void Foam::functionObjects::zetaSampler::writeZetaPlanesFile() const
{
    if (!Pstream::master()) return;
    if (zetaPlanesFileWritten_) return;

    std::ofstream os((outDirBase_ / "zetaPlanes.txt").c_str());
    os << "#==============================================================\n"
       << "# zetaSampler output spec\n"
       << "#==============================================================\n"
       << "#\n"
       << "# Per output time, this functionObject writes a directory:\n"
       << "#   postProcessing/" << "<this fO name>" << "/<time>/\n"
       << "# containing one binary file per MPI rank that produced any\n"
       << "# samples:\n"
       << "#   rank_<NNNNN>.bin                (NNNNN = zero-padded proc id)\n"
       << "#\n"
       << "# Each rank file is a stream of records (no per-file header):\n"
       << "#\n"
       << "#   record (variable size, repeated until EOF):\n"
       << "#     int64  levelIdx    (0 .. nZeta-1)\n"
       << "#     int64  fieldIdx    (index into 'fields' list below)\n"
       << "#     int64  i           (streamwise index, 0 .. Nxy(levelIdx)-1)\n"
       << "#     int64  j           (spanwise index,   0 .. Nxy(levelIdx)-1)\n"
       << "#     double val[ncomp(fieldIdx)]\n"
       << "#\n"
       << "# Nxy varies per level (see table below). ncomp varies per field\n"
       << "# (1=scalar, 3=vector, 9=tensor); see the 'fields' line. Postprocess\n"
       << "# reads each rank file, dispatches records by fieldIdx to the right\n"
       << "# (Nxy, Nxy, ncomp) NaN-filled array indexed [j, i, :].\n"
       << "#\n"
       << "# Cell centres on level k:\n"
       << "#   x_i = xmin + (i + 0.5) * (Lx / Nxy_k),  i = 0 .. Nxy_k-1\n"
       << "#   y_j = ymin + (j + 0.5) * (Ly / Nxy_k),  j = 0 .. Nxy_k-1\n"
       << "#\n"
       << "# Tensor component order (for val[c] when ncomp=9):\n"
       << "#   Txx, Txy, Txz, Tyx, Tyy, Tyz, Tzx, Tzy, Tzz\n"
       << "#\n"
       << "# Missing (i, j) pairs (target outside global domain or eta\n"
       << "# lookup failed) simply produce no record; the array entry stays\n"
       << "# NaN in postprocess.\n"
       << "#\n"
       << "# zeta-coordinate definition (Hara & Sullivan 2015 damped\n"
       << "# mapping):\n"
       << "#   z(x, y, zeta, t) = zeta + eta(x, y, t) * exp(-kp*|zeta|)\n"
       << "# Plane k samples z* = zeta_k + D_k * eta(x_i, y_j, t) with\n"
       << "# D_k = exp(-kp*|zeta_k|) (decay column in the table below),\n"
       << "# so planes follow the surface fully at zeta = 0 and relax to\n"
       << "# flat fixed-z planes over the scale 1/kp. kp = 0 recovers the\n"
       << "# uniform shift z* = eta + zeta_k. Metric terms for the\n"
       << "# postprocess (D = exp(-kp*|zeta|)):\n"
       << "#   dz/dzeta      = 1 - kp*sgn(zeta)*eta*D\n"
       << "#   dz/dx |_zeta  = D * deta/dx\n"
       << "#   dz/dt |_zeta  = D * deta/dt\n"
       << "#\n"
       << "# Per output time the master additionally writes the resampled\n"
       << "# surface elevation used to build the targets, one file per\n"
       << "# distinct tier resolution (wave-following mode only):\n"
       << "#   eta_<Nxy>.bin   (Nxy*Nxy doubles, row-major [j*Nxy + i],\n"
       << "#                    NaN where the eta lookup failed)\n"
       << "# These are the exact eta arrays the sampler used, so plane\n"
       << "# geometry and metric terms are reconstructible bit-for-bit\n"
       << "# offline without replicating the resampling convention.\n"
       << "#--------------------------------------------------------------\n"
       << "\n"
       << "xmin   " << std::scientific << std::setprecision(12) << xmin_ << "\n"
       << "ymin   " << ymin_ << "\n"
       << "Lx     " << Lx_   << "\n"
       << "Ly     " << Ly_   << "\n"
       << "kp     " << kp_   << "\n"
       << "nZeta  " << zetaLevels_.size() << "\n"
       << "fields";
    forAll(fieldNames_, i) os << " " << fieldNames_[i];
    os << "\n\n"
       << "# Per-level table:\n"
       << "# levelIndex   zetaValue [m]   Nxy   decay\n";
    forAll(zetaLevels_, k)
    {
        os << "  " << std::setw(4) << k
           << "    " << std::setw(20) << zetaLevels_[k]
           << "    " << std::setw(6)  << levelNxy_[k]
           << "    " << std::setw(20) << levelDecay_[k]
           << "\n";
    }
    os.close();

    zetaPlanesFileWritten_ = true;
}



// Walk from currentCelli along face neighbours toward target point pt,
// stopping when the cell containing pt is found or the walk exits the
// rank's local mesh. Returns -1 on failure, in which case caller should
// fall back to a global findCell.
//
// On a structured Cartesian-like AMR mesh with a vertical column of
// targets, this walk is at most a few face-hops per zeta increment.
// Sanitize a field name for use in a filesystem path. OpenFOAM derived
// fields like grad(U), grad(p) carry parentheses, dots, colons that some
// tools choke on. We replace any of those with '_'.
static Foam::word sanitizeForPath(const Foam::word& name)
{
    Foam::string s(name);
    for (size_t i = 0; i < s.size(); ++i)
    {
        const char c = s[i];
        if (c == '(' || c == ')' || c == '.' || c == ':' ||
            c == '/' || c == '\\' || c == ' ' || c == '\t')
        {
            s[i] = '_';
        }
    }
    return Foam::word(s, false);
}


static Foam::label advanceCellTowards
(
    const Foam::polyMesh& mesh,
    Foam::label currentCelli,
    const Foam::point& pt,
    Foam::label maxHops = 32
)
{
    if (currentCelli < 0 || currentCelli >= mesh.nCells()) return -1;

    const Foam::labelUList& own = mesh.faceOwner();
    const Foam::labelUList& nei = mesh.faceNeighbour();
    const Foam::cellList& cells = mesh.cells();
    const Foam::vectorField& Cf = mesh.faceCentres();
    const Foam::vectorField& Sf = mesh.faceAreas();

    Foam::label celli = currentCelli;
    for (Foam::label hop = 0; hop < maxHops; ++hop)
    {
        // Check whether pt is on the same side as cellCentre for every
        // face of the cell. If yes, pt is inside this cell.
        const Foam::cell& c = cells[celli];
        bool inside = true;
        Foam::label exitFace = -1;
        Foam::scalar mostExitDot = 0;

        forAll(c, fi)
        {
            const Foam::label faceI = c[fi];
            // Outward normal: if own[faceI] == celli, Sf points outward.
            // If nei[faceI] == celli (internal face), Sf points inward.
            Foam::vector nOut = Sf[faceI];
            if (faceI < mesh.nInternalFaces() && nei[faceI] == celli)
                nOut = -nOut;

            const Foam::scalar d = (pt - Cf[faceI]) & nOut;
            if (d > 0)
            {
                inside = false;
                if (d > mostExitDot)
                {
                    mostExitDot = d;
                    exitFace = faceI;
                }
            }
        }

        if (inside) return celli;
        if (exitFace < 0) return -1;
        if (exitFace >= mesh.nInternalFaces()) return -1;  // hit boundary

        // Cross to the neighbour through exitFace
        const Foam::label nbr =
            (own[exitFace] == celli) ? nei[exitFace] : own[exitFace];
        if (nbr < 0) return -1;
        celli = nbr;
    }
    return -1;
}


// Single-level streaming sampler. Caller supplies:
//  - per-level (nxy x nxy) eta lookup (already on this nxy grid)
//  - colCell cache from the *previous* level on the same nxy grid;
//    must be size nxy*nxy and pre-initialised to -1 the first time
//    the caller uses a new nxy.
// Master appends the resulting plane (size nxy*nxy*ncomp doubles)
// to os.
// Single-level streaming sampler with pre-built interpolator. The
// interpolator is constructed ONCE per field per execute() outside
// this template, then passed in unchanged across all zeta levels.
// This avoids rebuilding the cell-point interpolation table for every
// (level, field) pair, which was the dominant cost in the original
// per-level construction.
// Globals to accumulate per-execute time breakdown
static double       g_tFindHits  = 0;   // walk-from-seed succeeded
static double       g_tFindMiss  = 0;   // findCell fallback fired
static Foam::label  g_nFindMiss  = 0;
static Foam::label  g_nFindHits  = 0;
static double       g_tInterp    = 0;
static double       g_tReduce    = 0;
static double       g_tWrite     = 0;

// Per-rank sample-and-write. Each rank writes its own hits inline to
// its own file. No gather, no master serialization, no per-(level,field)
// collective. Record layout (one per hit):
//   int64 levelIdx, int64 fieldIdx, int64 i, int64 j
//   double val[ncomp]
// ncomp varies by field (1, 3, or 9), so records are variable length;
// the postprocess reader uses fieldIdx to look up ncomp from manifest.txt.
//
// The target height at level k follows the Hara-Sullivan damped
// mapping: z* = zeta_k + decay_k * eta_ij with decay_k =
// exp(-kp*|zeta_k|) precomputed by the caller. decay_k = 1 (kp = 0)
// recovers the uniform shift z* = eta_ij + zeta_k.
template<class Type>
static void samplePerLevelToRankFile
(
    const Foam::fvMesh& mesh,
    const Foam::interpolationCellPoint<Type>& interp,
    const Foam::List<Foam::scalar>& etaLocal,   // size nxy*nxy
    const Foam::List<Foam::scalar>& zRange,     // size 2*nxy*nxy: [zMin,zMax]
    Foam::label nxy,
    Foam::scalar zeta_k,
    Foam::scalar decay_k,                       // D_k = exp(-kp*|zeta_k|)
    Foam::scalar xmin, Foam::scalar ymin,
    Foam::scalar dx,   Foam::scalar dy,
    Foam::labelList& colCell,                   // size nxy*nxy
    const Foam::boundBox& rankBBpad,
    Foam::label levelIdx,
    Foam::label fieldIdx,
    Foam::label& nHitsThisLF,                   // hit counter for this (level,field)
    std::ofstream& rankOs                       // append-mode rank file
)
{
    using namespace Foam;

    const label ncomp = pTraits<Type>::nComponents;

    for (label j = 0; j < nxy; ++j)
    for (label i = 0; i < nxy; ++i)
    {
        const scalar x = xmin + (i + 0.5) * dx;
        const scalar y = ymin + (j + 0.5) * dy;
        const scalar eta_ij = etaLocal[j * nxy + i];
        if (std::isnan(eta_ij)) continue;

        const scalar zStar = eta_ij * decay_k + zeta_k;
        const point pt(x, y, zStar);

        if (!rankBBpad.contains(pt)) continue;

        // z-range filter: this rank's cells at (i, j) span [zMin, zMax].
        // If z* is outside that range no cell can contain it on this
        // rank. Removes ~95% of findCell misses caused by other-rank
        // ownership in z. SMALL pad covers floating-point edges.
        const label b = j * nxy + i;
        const scalar zMinIj = zRange[2*b];
        const scalar zMaxIj = zRange[2*b + 1];
        if (zStar < zMinIj - SMALL || zStar > zMaxIj + SMALL) continue;

        // Walk from cached column cell, fall back to findCell.
        label celli = -1;
        const label seed = colCell[j * nxy + i];
        const double t0a = mesh.time().elapsedCpuTime();
        if (seed >= 0)
        {
            celli = advanceCellTowards(mesh, seed, pt);
        }
        bool walkHit = (celli >= 0);
        if (celli < 0)
        {
            celli = mesh.findCell(pt);
        }
        const double t0b = mesh.time().elapsedCpuTime();
        if (walkHit) { g_tFindHits += (t0b - t0a); g_nFindHits++; }
        else         { g_tFindMiss += (t0b - t0a); g_nFindMiss++; }

        colCell[j * nxy + i] = celli;

        if (celli < 0) continue;

        const double t1a = mesh.time().elapsedCpuTime();
        const Type val = interp.interpolate(pt, celli);
        const double t1b = mesh.time().elapsedCpuTime();
        g_tInterp += (t1b - t1a);

        // Write record directly to this rank's file.
        const double t3a = mesh.time().elapsedCpuTime();
        const int64_t hdr[4] =
        {
            int64_t(levelIdx),
            int64_t(fieldIdx),
            int64_t(i),
            int64_t(j)
        };
        rankOs.write
        (
            reinterpret_cast<const char*>(hdr),
            std::streamsize(sizeof(hdr))
        );
        scalar buf[9];   // max ncomp is 9 (tensor)
        for (label c = 0; c < ncomp; ++c)
        {
            buf[c] = component(val, direction(c));
        }
        rankOs.write
        (
            reinterpret_cast<const char*>(buf),
            std::streamsize(ncomp * sizeof(scalar))
        );
        const double t3b = mesh.time().elapsedCpuTime();
        g_tWrite += (t3b - t3a);

        nHitsThisLF++;
    }
}


void Foam::functionObjects::zetaSampler::resampleEtaToLevel
(
    const List<scalar>& etaSrc,
    label NxSrc, label NySrc,
    scalar dlSrc, scalar xSrc0, scalar ySrc0,
    label nxy,
    List<scalar>& etaOut
) const
{
    const scalar nanVal = std::numeric_limits<scalar>::quiet_NaN();
    etaOut.setSize(nxy * nxy, nanVal);

    const scalar dx = Lx_ / nxy;
    const scalar dy = Ly_ / nxy;

    for (label j = 0; j < nxy; ++j)
    for (label i = 0; i < nxy; ++i)
    {
        const scalar x = xmin_ + (i + 0.5) * dx;
        const scalar y = ymin_ + (j + 0.5) * dy;
        const label is = label(round((x - xSrc0) / dlSrc));
        const label js = label(round((y - ySrc0) / dlSrc));
        if (is < 0 || is >= NxSrc || js < 0 || js >= NySrc) continue;
        etaOut[j * nxy + i] = etaSrc[js * NxSrc + is];
    }
}


bool Foam::functionObjects::zetaSampler::execute()
{
    Info<< "zetaSampler::execute()  t=" << mesh_.time().timeName()
        << "  fieldsChecked=" << fieldsChecked_
        << endl;

    // First-call field-existence check. Aborts if any field in the
    // 'fields' list is not registered. Deferred from read() because
    // some producing functionObjects (e.g. grad) only register their
    // output field on their own first execute(), not in their
    // constructor. By the time we reach this point in execute(), all
    // earlier-listed functionObjects have run their first execute(),
    // so the fields they produce should now be in the registry.
    if (!fieldsChecked_)
    {
        wordList missing;
        forAll(fieldNames_, fi)
        {
            const word& fname = fieldNames_[fi];
            const bool found =
                   mesh_.foundObject<volScalarField>(fname)
                || mesh_.foundObject<volVectorField>(fname)
                || mesh_.foundObject<volTensorField>(fname);
            if (!found) missing.append(fname);
        }

        if (!missing.empty())
        {
            // Fail HARD - bypass functionObject loader catch by calling
            // Foam::FatalError + abort(). FatalErrorInFunction with
            // exit(FatalError) inside a functionObject is intercepted
            // by the runtime selection wrapper and downgraded to a
            // warning that lets the solver continue. We don't want
            // that: missing fields means meaningless output.
            FatalErrorInFunction
                << nl
                << "================================================" << nl
                << "zetaSampler ABORT: missing fields in registry"     << nl
                << "================================================" << nl
                << "  Missing: " << missing << nl
                << "  Action:  add a producing functionObject earlier" << nl
                << "           in controlDict.functions, with the same"<< nl
                << "           executeInterval as zetaSampler. For"    << nl
                << "           example, 'grad' fO produces grad(<F>)." << nl
                << "  Note:    these checks fire on the first execute,"<< nl
                << "           not at construction, because OpenFOAM"  << nl
                << "           derived-field producers only register"  << nl
                << "           their output on first execute()."       << nl
                << "================================================" << nl
                << exit(FatalError);

            // Belt-and-braces: if the loader catches the FatalError
            // above, kill the solver outright.
            Foam::FatalError.exit(1);
        }
        fieldsChecked_ = true;
    }

    // Look up the binned eta grid (IOList<scalar>) and its spec
    // (IOdictionary) from the mesh objectRegistry. Both are registered
    // by the surfaceElevation functionObject under fixed names. If
    // either is missing, surfaceElevation has not run yet (or the
    // user named it something other than etaSourceName_); abort.
    const word etaListName = etaSourceName_ + ":etaGrid";
    const word etaSpecName = etaSourceName_ + ":etaGridSpec";

    if
    (
        !mesh_.foundObject<IOList<scalar>>(etaListName)
     || !mesh_.foundObject<IOdictionary>(etaSpecName)
    )
    {
        Info<< "zetaSampler: eta registry entries not yet present "
            << "(looking for '" << etaListName << "' and '"
            << etaSpecName << "'); skipping this step." << endl;
        return true;
    }

    const IOList<scalar>&  etaGridRefIO =
        mesh_.lookupObject<IOList<scalar>>(etaListName);
    const IOdictionary&    etaSpec =
        mesh_.lookupObject<IOdictionary>(etaSpecName);

    const List<scalar>& etaGridRef = etaGridRefIO;
    const label  srcNx = etaSpec.get<label>("Nx");
    const label  srcNy = etaSpec.get<label>("Ny");
    const scalar srcDl = etaSpec.get<scalar>("dl");
    const scalar srcX0 = etaSpec.get<scalar>("xmin");
    const scalar srcY0 = etaSpec.get<scalar>("ymin");

    // Monotonicity check of the damped mapping on the instantaneous
    // eta grid. The mapping z = zeta + eta*exp(-kp*|zeta|) folds when
    // kp*|eta| >= 1: the Jacobian dz/dzeta = 1 - kp*sgn(zeta)*eta*D
    // reaches its minimum 1 - kp*max|eta| at the surface (over crests
    // on the air side, under troughs on the water side). A folded
    // mapping inverts the vertical ordering of planes with no other
    // symptom, so this is fatal, with the same hard-abort pattern as
    // the field check. The minimum instantaneous Jacobian is reported
    // each execute so marginal snapshots are visible in the log.
    if (!flatInterface_ && kp_ > 0)
    {
        scalar maxAbsEta = 0;
        forAll(etaGridRef, n)
        {
            const scalar e = etaGridRef[n];
            if (!std::isnan(e))
            {
                const scalar a = mag(e);
                if (a > maxAbsEta) maxAbsEta = a;
            }
        }
        reduce(maxAbsEta, maxOp<scalar>());

        const scalar minJ = 1 - kp_ * maxAbsEta;

        Info<< "zetaSampler: mapping monotonicity"
            << "  max|eta|=" << maxAbsEta
            << "  kp*max|eta|=" << kp_ * maxAbsEta
            << "  min(dz/dzeta)=" << minJ << endl;

        if (minJ <= 0)
        {
            FatalErrorInFunction
                << nl
                << "================================================" << nl
                << "zetaSampler ABORT: damped mapping has folded"      << nl
                << "================================================" << nl
                << "  kp           = " << kp_ << nl
                << "  max|eta|     = " << maxAbsEta << nl
                << "  kp*max|eta|  = " << kp_ * maxAbsEta
                << "  (must be < 1)" << nl
                << "  A folded mapping places a nominally higher zeta" << nl
                << "  level at lower physical z over the steepest"     << nl
                << "  crest/trough, silently inverting the vertical"   << nl
                << "  ordering of the planes. Reduce kp."              << nl
                << "================================================" << nl
                << exit(FatalError);

            // Belt-and-braces: if the loader catches the FatalError
            // above, kill the solver outright.
            Foam::FatalError.exit(1);
        }
    }

    writeZetaPlanesFile();

    // Per-time output directory
    const fileName timeDir = outDirBase_ / mesh_.time().timeName();
    // Master creates the time directory; a sync (cheap scalar reduce)
    // ensures every rank sees it before opening its rank-local file.
    if (Pstream::master()) mkDir(timeDir);
    {
        scalar dummy = 0;
        reduce(dummy, sumOp<scalar>());
    }

    const scalar t0 = mesh_.time().elapsedCpuTime();

    g_tFindHits = g_tFindMiss = g_tInterp = g_tReduce = g_tWrite = 0;
    g_nFindHits = g_nFindMiss = 0;

    // Cache per-tier resampled eta and colCell. Keys are nxy values.
    // Pre-build for every distinct nxy used across all levels so the
    // outer field loop has no nxy bookkeeping inside. The per-level
    // decay factor D_k is NOT folded into these arrays: the cache is
    // keyed by nxy and shared across all levels in a tier, while D_k
    // differs per level. D_k is applied in the sampling loop instead.
    Map<List<scalar>>  etaPerNxy;
    Map<labelList>     colCellPerNxy;

    forAll(zetaLevels_, k)
    {
        const label nxy = levelNxy_[k];
        if (!etaPerNxy.found(nxy))
        {
            List<scalar> etaTmp;
            if (flatInterface_)
            {
                // Rigid-lid: target z* = 0 + zeta_k, so etaLocal is zero.
                etaTmp.setSize(nxy * nxy, 0.0);
            }
            else
            {
                resampleEtaToLevel
                (
                    etaGridRef, srcNx, srcNy, srcDl, srcX0, srcY0,
                    nxy, etaTmp
                );
            }
            etaPerNxy.insert(nxy, etaTmp);
            colCellPerNxy.insert(nxy, labelList(nxy * nxy, -1));
        }
    }

    // Write the resampled eta arrays the sampler actually uses, one
    // binary per distinct tier resolution per snapshot (master only,
    // wave-following mode only). Layout: nxy*nxy doubles, row-major
    // [j*nxy + i], NaN where the lookup failed. Under the damped
    // mapping the plane geometry and every metric term (dz/dzeta,
    // dz/dx, dz/dt) depend on eta scaled per level by D_k, so an
    // offline reconstruction of eta from the surfaceElevation output
    // would have to replicate this resampling bit-for-bit; any
    // mismatch enters the metrics scaled by D_k, varying across
    // levels. Writing the arrays used removes that failure class.
    if (!flatInterface_ && Pstream::master())
    {
        forAllConstIters(etaPerNxy, it)
        {
            const label nxy = it.key();
            const List<scalar>& etaTmp = it.val();

            std::ostringstream ess;
            ess << "eta_" << nxy << ".bin";
            const fileName etaFile = timeDir / ess.str();

            std::ofstream eos(etaFile.c_str(), std::ios::binary);
            if (!eos.is_open())
            {
                FatalErrorInFunction
                    << "zetaSampler: failed to open eta output file "
                    << etaFile << exit(FatalError);
            }
            eos.write
            (
                reinterpret_cast<const char*>(etaTmp.cdata()),
                std::streamsize(etaTmp.size() * sizeof(scalar))
            );
            eos.close();
        }
    }

    // Per-tier 2D z-range map of rank-local cell centres. zRangePerNxy
    // [nxy] is a List<scalar> of size 2*nxy*nxy holding [zMin_ij, zMax_ij]
    // pairs in row-major order (j*nxy + i). For a sample at (i, j) on
    // tier nxy with target z*, if z* is outside this rank's local
    // [zMin_ij, zMax_ij] at that bin, no cell on this rank can contain
    // it - skip findCell entirely. This is the z-band generalization
    // of the rank AABB and removes ~95% of the findCell misses.
    Map<List<scalar>> zRangePerNxy;
    {
        const vectorField& C = mesh_.C().internalField();
        const cellList& cells = mesh_.cells();
        const faceList& faces = mesh_.faces();
        const pointField& pts = mesh_.points();
        forAllConstIters(etaPerNxy, it)
        {
            const label nxy = it.key();
            const scalar dx = Lx_/nxy;
            const scalar dy = Ly_/nxy;

            const scalar zHi =  std::numeric_limits<scalar>::infinity();
            const scalar zLo = -std::numeric_limits<scalar>::infinity();
            List<scalar> zr(2 * nxy * nxy);
            for (label b = 0; b < nxy * nxy; ++b)
            {
                zr[2*b]   =  zHi;   // zMin starts at +inf
                zr[2*b+1] =  zLo;   // zMax starts at -inf
            }

            for (label celli = 0; celli < C.size(); ++celli)
            {
                // Walk the cell's vertices to get its full (x,y,z) extent.
                // Register the cell's z-range at every (i,j) bin its
                // (x,y) footprint touches, not just the centre's bin.
                // AMR-fine cells (1.25 mm) on a 1.25 mm sampler tier can
                // straddle bin boundaries; missing those would leave
                // valid samples unreachable on that rank.
                scalar xCellMin =  zHi, xCellMax =  zLo;
                scalar yCellMin =  zHi, yCellMax =  zLo;
                scalar zCellMin =  zHi, zCellMax =  zLo;
                const labelList& cFaces = cells[celli];
                forAll(cFaces, fi)
                {
                    const face& f = faces[cFaces[fi]];
                    forAll(f, pi)
                    {
                        const point& p = pts[f[pi]];
                        if (p.x() < xCellMin) xCellMin = p.x();
                        if (p.x() > xCellMax) xCellMax = p.x();
                        if (p.y() < yCellMin) yCellMin = p.y();
                        if (p.y() > yCellMax) yCellMax = p.y();
                        if (p.z() < zCellMin) zCellMin = p.z();
                        if (p.z() > zCellMax) zCellMax = p.z();
                    }
                }

                label iLo = label(std::floor((xCellMin - xmin_) / dx));
                label iHi = label(std::floor((xCellMax - xmin_) / dx));
                label jLo = label(std::floor((yCellMin - ymin_) / dy));
                label jHi = label(std::floor((yCellMax - ymin_) / dy));
                if (iLo < 0)        iLo = 0;
                if (iHi >= nxy)     iHi = nxy - 1;
                if (jLo < 0)        jLo = 0;
                if (jHi >= nxy)     jHi = nxy - 1;
                if (iLo > iHi || jLo > jHi) continue;

                for (label j = jLo; j <= jHi; ++j)
                for (label i = iLo; i <= iHi; ++i)
                {
                    const label b = j * nxy + i;
                    if (zCellMin < zr[2*b])   zr[2*b]   = zCellMin;
                    if (zCellMax > zr[2*b+1]) zr[2*b+1] = zCellMax;
                }
            }
            zRangePerNxy.insert(nxy, zr);
        }
    }

    // Outer loop: zeta levels. The colCell cache persists between
    // adjacent levels within a tier so the cell-walk seeds neighbour
    // traversal (each level's findCell results seed the next level).
    // Inner loop: fields. Each (level, field) constructs its own
    // interpolator inline; this is the structure that worked before
    // we tried (and failed) to invert the loops to share interpolators.
    Info<< "zetaSampler: walking " << zetaLevels_.size()
        << " levels x " << fieldNames_.size() << " fields" << endl;

    // Rank-local AABB, computed ONCE per execute. mesh.bounds() returns
    // the GLOBAL box (already reduced across procs), which would skip
    // nothing. Build the local-only AABB from this rank's points.
    boundBox rankBBpad;
    {
        const pointField& lp = mesh_.points();
        if (lp.empty())
        {
            rankBBpad = boundBox(point::zero, point::zero);
        }
        else
        {
            point lo(lp[0]), hi(lp[0]);
            for (label p = 1; p < lp.size(); ++p)
            {
                lo = min(lo, lp[p]);
                hi = max(hi, lp[p]);
            }
            const vector pad(SMALL, SMALL, SMALL);
            rankBBpad = boundBox(lo - pad, hi + pad);
        }
    }

    // Build interpolators ONCE before the level loop. Each
    // interpolationCellPoint<Type> constructor builds the cell-to-point
    // table for the entire field; doing this per (level, field) is
    // ~636 builds per execute. With prebuild, it's just N_fields.
    enum FieldType { TYPE_SCALAR, TYPE_VECTOR, TYPE_TENSOR, TYPE_NONE };
    List<FieldType> ftype(fieldNames_.size(), TYPE_NONE);
    PtrList<interpolationCellPoint<scalar>> interpS(fieldNames_.size());
    PtrList<interpolationCellPoint<vector>> interpV(fieldNames_.size());
    PtrList<interpolationCellPoint<tensor>> interpT(fieldNames_.size());

    forAll(fieldNames_, fi)
    {
        const word& fname = fieldNames_[fi];
        if (mesh_.foundObject<volScalarField>(fname))
        {
            ftype[fi] = TYPE_SCALAR;
            interpS.set
            (
                fi,
                new interpolationCellPoint<scalar>
                (
                    mesh_.lookupObject<volScalarField>(fname)
                )
            );
        }
        else if (mesh_.foundObject<volVectorField>(fname))
        {
            ftype[fi] = TYPE_VECTOR;
            interpV.set
            (
                fi,
                new interpolationCellPoint<vector>
                (
                    mesh_.lookupObject<volVectorField>(fname)
                )
            );
        }
        else if (mesh_.foundObject<volTensorField>(fname))
        {
            ftype[fi] = TYPE_TENSOR;
            interpT.set
            (
                fi,
                new interpolationCellPoint<tensor>
                (
                    mesh_.lookupObject<volTensorField>(fname)
                )
            );
        }
    }

    // Each rank opens ONE binary file and writes its hits across all
    // levels and all fields inline (no gather, no master serialization).
    // Each record carries (levelIdx, fieldIdx, i, j, val[ncomp(fieldIdx)]).
    // Postprocess reads all rank files and reassembles per-level
    // per-field arrays using ncomp lookup from manifest.txt.
    std::ofstream rankOs;
    {
        std::ostringstream rss;
        rss << "rank_" << std::setfill('0') << std::setw(5)
            << Pstream::myProcNo() << ".bin";
        const fileName rankFile = timeDir / rss.str();
        rankOs.open(rankFile.c_str(), std::ios::binary);
        if (!rankOs.is_open())
        {
            FatalErrorInFunction
                << "zetaSampler: failed to open per-rank output file "
                << rankFile << exit(FatalError);
        }
    }

    // Per-(level, field) hit counter on each rank. Reduced across ranks
    // at the end of execute() so the master can report NaN coverage
    // without needing a postprocessor. nHits[k * nF + fi] = total hits
    // at level k for field fi across all ranks.
    const label nL = zetaLevels_.size();
    const label nF = fieldNames_.size();
    labelList nHitsPerLF(nL * nF, 0);

    label nLevelsDone = 0;
    scalar tBlockStart = mesh_.time().elapsedCpuTime();
    forAll(zetaLevels_, k)
    {
        const label  nxy     = levelNxy_[k];
        const scalar zeta_k  = zetaLevels_[k];
        const scalar decay_k = levelDecay_[k];
        const List<scalar>& etaLocal = etaPerNxy[nxy];
        const List<scalar>& zRange   = zRangePerNxy[nxy];
        labelList&          colCell  = colCellPerNxy[nxy];

        forAll(fieldNames_, fi)
        {
            switch (ftype[fi])
            {
                case TYPE_SCALAR:
                    samplePerLevelToRankFile<scalar>
                    (
                        mesh_, interpS[fi], etaLocal, zRange,
                        nxy, zeta_k, decay_k,
                        xmin_, ymin_, Lx_/nxy, Ly_/nxy,
                        colCell, rankBBpad, k, fi, nHitsPerLF[k * nF + fi], rankOs
                    );
                    break;
                case TYPE_VECTOR:
                    samplePerLevelToRankFile<vector>
                    (
                        mesh_, interpV[fi], etaLocal, zRange,
                        nxy, zeta_k, decay_k,
                        xmin_, ymin_, Lx_/nxy, Ly_/nxy,
                        colCell, rankBBpad, k, fi, nHitsPerLF[k * nF + fi], rankOs
                    );
                    break;
                case TYPE_TENSOR:
                    samplePerLevelToRankFile<tensor>
                    (
                        mesh_, interpT[fi], etaLocal, zRange,
                        nxy, zeta_k, decay_k,
                        xmin_, ymin_, Lx_/nxy, Ly_/nxy,
                        colCell, rankBBpad, k, fi, nHitsPerLF[k * nF + fi], rankOs
                    );
                    break;
                case TYPE_NONE:
                    break;
            }
        }

        nLevelsDone++;
        if (nLevelsDone % 10 == 0
         || nLevelsDone == zetaLevels_.size())
        {
            scalar tNow = mesh_.time().elapsedCpuTime();
            scalar dtBlock = tNow - tBlockStart;
            scalar dtTotal = tNow - t0;
            tBlockStart = tNow;
            reduce(dtBlock, maxOp<scalar>());
            reduce(dtTotal, maxOp<scalar>());
            Info<< "zetaSampler:  level " << nLevelsDone
                << "/" << zetaLevels_.size()
                << "  zeta=" << zeta_k
                << "  nxy=" << nxy
                << "  dt(last 10)=" << dtBlock << "s"
                << "  total=" << dtTotal << "s" << endl;
        }
    }

    rankOs.close();


    const scalar t1 = mesh_.time().elapsedCpuTime();
    scalar dt = t1 - t0;
    reduce(dt, maxOp<scalar>());

    Info<< "zetaSampler: t=" << mesh_.time().timeName()
        << "  T(sample+IO)=" << dt << "s" << endl;

    // Per-bucket timing breakdown (rank-local; reduce to master)
    {
        scalar tFindHits = g_tFindHits, tFindMiss = g_tFindMiss;
        scalar tInterp   = g_tInterp,   tReduce   = g_tReduce;
        scalar tWrite    = g_tWrite;
        label  nFindHits = g_nFindHits, nFindMiss = g_nFindMiss;

        reduce(tFindHits, sumOp<scalar>());
        reduce(tFindMiss, sumOp<scalar>());
        reduce(tInterp,   sumOp<scalar>());
        reduce(tReduce,   sumOp<scalar>());
        reduce(tWrite,    sumOp<scalar>());
        reduce(nFindHits, sumOp<label>());
        reduce(nFindMiss, sumOp<label>());

        Info<< "zetaSampler: profile (sum across ranks)" << nl
            << "  T(walk hits) =" << tFindHits
            << " s   nHits=" << nFindHits << nl
            << "  T(findCell)  =" << tFindMiss
            << " s   nMiss=" << nFindMiss << nl
            << "  T(interp)    =" << tInterp << " s" << nl
            << "  T(reduce)    =" << tReduce << " s" << nl
            << "  T(write)     =" << tWrite  << " s"
            << endl;
    }

    // Reduce per-(level, field) hit counts to master and report.
    Pstream::listReduce(nHitsPerLF, sumOp<label>());
    if (Pstream::master())
    {
        label totalExpected = 0;
        label totalHits     = 0;
        label nUnder        = 0;   // planes with missing samples (NaN)
        label nOver         = 0;   // planes with duplicate samples
        for (label k = 0; k < nL; ++k)
        {
            const label nxy = levelNxy_[k];
            const label exp1 = nxy * nxy;
            for (label fi = 0; fi < nF; ++fi)
            {
                const label hits = nHitsPerLF[k * nF + fi];
                totalExpected += exp1;
                totalHits     += hits;
                if (hits < exp1) nUnder++;
                if (hits > exp1) nOver++;
            }
        }

        Info<< "zetaSampler: sample coverage report" << nl;
        Info<< "  expected samples (all planes): " << totalExpected << nl;
        Info<< "  recorded samples:              " << totalHits << nl;

        if (totalHits == totalExpected && nUnder == 0 && nOver == 0)
        {
            Info<< "  -> every plane fully covered, no duplicates." << endl;
        }
        else
        {
            if (nUnder > 0)
            {
                Info<< "  WARNING: " << nUnder << " of " << (nL*nF)
                    << " planes have MISSING samples (NaN holes)" << endl;
            }
            if (nOver > 0)
            {
                Info<< "  WARNING: " << nOver << " of " << (nL*nF)
                    << " planes have DUPLICATE samples (same (i,j) "
                    << "written by more than one rank)" << endl;
            }
        }
    }

    return true;
}


bool Foam::functionObjects::zetaSampler::write()
{
    return true;
}
