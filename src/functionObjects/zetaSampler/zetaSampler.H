/*---------------------------------------------------------------------------*\
            Copyright (c) 2026, Asim Onder
-------------------------------------------------------------------------------
License
    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY. See the GNU General Public License for more details.

Class
    Foam::functionObjects::zetaSampler

Description
    Samples a list of registered volume fields onto a regular Cartesian
    (x_i, y_j) grid at user-specified zeta levels above and below the
    instantaneous wave surface eta(x_i, y_j, t). Per output time, writes
    one binary file per rank with the rank's contributions to the global
    Nx*Ny*Nzeta grid; the master also writes a header file with the
    sampling spec and the per-tier resampled eta arrays.

    The wave-following coordinate uses the exponentially damped mapping
    of Hara & Sullivan (2015),

        z(x, y, zeta, t) = zeta + eta(x, y, t) * exp(-kp * |zeta|)

    so that planes follow the surface fully at zeta = 0 and relax to
    flat fixed-z planes over the vertical scale 1/kp. The target
    sampling point at horizontal index (i, j) and zeta level zeta_k is

        z*(i, j, k, t) = zeta_k + D_k * eta(x_i, y_j, t),
        D_k = exp(-kp * |zeta_k|)

    with the per-level decay factors D_k precomputed at construction.
    Setting kp = 0 (the default) gives D_k = 1 for all levels and
    recovers the uniform-shift mapping z* = eta + zeta_k exactly.

    The mapping is monotone in zeta provided kp * max|eta| < 1. This is
    checked on the instantaneous eta grid at every execute and the
    solver is aborted on violation, because a folded mapping silently
    inverts the vertical ordering of the planes. The same factor is the
    Jacobian dz/dzeta = 1 - kp * sgn(zeta) * eta * D needed by the
    postprocess metric corrections; the minimum instantaneous Jacobian
    is reported each execute.

    For each target, OpenFOAM's interpolationCellPoint locates the AMR
    cell containing the point and returns a trilinear interpolant from
    cell-center values. AMR refinement transitions are handled natively
    by the polyMesh.

    The eta(x, y, t) field is consumed from the surfaceElevation
    functionObject's broadcast grid. Both functionObjects must be active
    in the same controlDict; zetaSampler must run at an executeInterval
    that is a multiple of surfaceElevation's executeInterval, so the
    eta field is current.

Usage
    Add a `zetaSampler` block to `system/controlDict.functions`. The
    `surfaceElevation` functionObject must also be active and run at an
    `executeInterval` that divides this object's `executeInterval`, so a
    fresh eta(x,y,t) is available on every call. Any field listed in
    `fields` must be registered in the mesh objectRegistry by the time
    this functionObject is constructed; produce derived fields (e.g.
    grad(U), grad(p)) via earlier functionObjects in the same controlDict.

    Example:
    \verbatim
    zetaSampler
    {
        type            zetaSampler;
        libs            ("libgeometricVoFCartesian.so");
        executeControl  adjustableRunTime;
        executeInterval 0.046;                  // T_p/8 for Case 1

        // Source of eta(x,y,t) (must run earlier in controlDict and at
        // an executeInterval that divides this object's executeInterval).
        etaSource       surfaceElevation;

        // Decay rate of the wave-following deformation (Hara &
        // Sullivan 2015), typically the spectral peak wavenumber.
        // Planes follow the surface fully at zeta = 0 and relax to
        // flat fixed-z planes over the scale 1/kp. Monotonicity
        // requires kp * max|eta| < 1; checked fatally each execute.
        // Optional; default 0 recovers the uniform shift
        // z* = eta + zeta_k.
        kp              32.0;                   // rad/m, Case 1

        // Optional. If true, ignore eta entirely - sample on flat
        // horizontal planes z = zeta_k instead of wave-following
        // z = zeta_k + eta * exp(-kp |zeta_k|). Use for rigid-lid
        // horizontal-plane statistics (Favre averaging,
        // oceanographic-style closures). Default: false.
        flatInterface   false;

        // Zeta level distribution. Built internally as a geometric
        // series outward from zeta = 0 on each side. The first non-zero
        // level on each side sits at +/- dz0; subsequent spacings are
        // multiplied by stretch. Walking stops on each side when the
        // next |zeta| would exceed zMaxAir / zMaxWater. Asymmetric
        // extents are natural (air-side wave layer is typically thinner
        // than water-side WBL + drift region).
        dz0             0.0004;                 // 0.4 mm at surface
        stretch         1.10;                   // 10% per level
        zMaxAir         0.040;                  // ~1.2/k_p, Case 1
        zMaxWater       0.160;                  // ~4.7/k_p, deeper WBL

        // Horizontal extent: taken automatically from mesh.bounds().
        // No dict entries needed.

        // Tier lists: separate lists for air side (zeta > 0) and water
        // side (zeta < 0). Each entry is (Nxy  |zeta|_max). A tier
        // spans from the previous tier's |zeta|_max (or 0 for the
        // innermost) up to its own. Order in the dict doesn't matter;
        // each list is sorted ascending by |zeta|_max at read time.
        // The outermost tier of each list must reach at least zMaxAir
        // and zMaxWater respectively, or sampler aborts at startup.
        tiersAir
        (
            (1024  0.012 )                      // air-side surface tier
            ( 512  0.040 )
            ( 256  0.060 )                      // air-side outer tier
        );
        tiersWater
        (
            (1024  0.012 )                      // water-side surface tier
            ( 512  0.060 )
            ( 256  0.160 )                      // water-side outer tier
        );

        // Fields to sample. Each must be registered in the mesh
        // objectRegistry as a vol(Scalar|Vector|Tensor)Field. Type
        // dispatch is automatic. Construction aborts if any name is
        // not found - no warnings, no silent skipping.
        fields          (U gradU p);
    }
    \endverbatim

Output
    postProcessing/<name>/
        zetaPlanes.txt                          // spec (written once)
        zetaLevels.txt                          // levels + Nxy + D_k
        <time1>/
            rank_00000.bin                      // per-rank sample records
            rank_00001.bin
            ...
            eta_1024.bin                        // resampled eta per tier
            eta_512.bin                         //   resolution (master,
            eta_256.bin                         //   wave-following mode)
        <time2>/
            ...

    Each rank file is a stream of variable-length records
    (levelIdx, fieldIdx, i, j, val[ncomp]); see zetaPlanes.txt for the
    full layout convention, the (levelIndex, zeta, Nxy, decay) table
    and the eta file layout. Each eta_<Nxy>.bin holds the exact
    Nxy*Nxy resampled surface elevation used to build the targets at
    that tier resolution (row-major [j*Nxy + i], NaN where the lookup
    failed), so plane geometry and metric terms are reconstructible
    bit-for-bit offline. zetaLevels.txt is a plain four-column text
    file (index, zeta, Nxy, decay) convenient for plotting the level
    distribution.

\*---------------------------------------------------------------------------*/

#ifndef zetaSampler_H
#define zetaSampler_H

#include "fvMeshFunctionObject.H"
#include "volFieldsFwd.H"
#include "scalarList.H"
#include "wordList.H"
#include "Tuple2.H"
#include <fstream>

namespace Foam
{

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                         Class zetaSampler Declaration
\*---------------------------------------------------------------------------*/

class zetaSampler
:
    public fvMeshFunctionObject
{
private:

    // Configuration

        //- Name of the surfaceElevation functionObject (eta source)
        word etaSourceName_;

        //- Stretching parameters for the zeta level distribution.
        //  Levels are built internally as a geometric series outward
        //  from zeta = 0 on each side independently:
        //    zeta_k+ = sum_{j=0..k-1} dz0 * stretch^j   (air side)
        //    zeta_k- = -sum_{j=0..k-1} dz0 * stretch^j  (water side)
        //  Stops on each side when |zeta| would exceed zMaxAir/zMaxWater.
        //  Asymmetric extents are natural since the air-side wave layer
        //  is typically thinner than the water-side WBL + drift region.
        scalar dz0_;
        scalar stretch_;
        scalar zMaxAir_;
        scalar zMaxWater_;

        //- Decay rate of the wave-following deformation (Hara &
        //  Sullivan 2015 damped mapping), typically the spectral peak
        //  wavenumber. The plane at level zeta_k deforms by
        //  D_k * eta(x, y, t) with D_k = exp(-kp * |zeta_k|).
        //  kp = 0 (default) recovers the uniform shift.
        scalar kp_;

        //- Horizontal extent (shared across all tiers)
        scalar xmin_;
        scalar ymin_;
        scalar Lx_;
        scalar Ly_;

        //- If true, ignore eta(x, y, t) entirely and treat the surface
        //  as a flat plane at z = 0. Sample targets become z* = zeta_k.
        //  Used for rigid-lid horizontal-plane statistics (Favre-style
        //  averaging across the full domain) as a complement to the
        //  wave-following mode. Default false.
        bool flatInterface_;

        //- Tier specification: separate lists for air side (zeta > 0)
        //  and water side (zeta < 0). Each entry is (N_xy, |zeta|_max).
        //  A tier spans from the previous tier's |zeta|_max (or 0 for
        //  the innermost) up to its own. Lookup matches a level to the
        //  first tier in the appropriate list whose |zeta|_max >=
        //  |zeta_k|. Lists are sorted ascending by |zeta|_max at read
        //  time so order in the dict doesn't matter. Asymmetry is
        //  natural: AMR refinement near the surface is typically
        //  asymmetric (e.g. extends further into the air than into
        //  the water, or vice versa).
        struct Tier
        {
            label  Nxy;
            scalar zMaxAbs;
        };
        List<Tier> tiersAir_;
        List<Tier> tiersWater_;

        //- Names of fields to sample (must exist in objectRegistry)
        wordList fieldNames_;


    // Cached state

        //- Sorted zeta levels (m)
        scalarList zetaLevels_;

        //- Per-level Nx (= Ny) assigned from tiers, sized to zetaLevels_
        labelList levelNxy_;

        //- Per-level decay factor D_k = exp(-kp * |zeta_k|), sized to
        //  zetaLevels_. Precomputed once after buildZetaLevels(); the
        //  target at level k is z* = zeta_k + D_k * eta(x_i, y_j, t).
        //  All ones when kp = 0.
        scalarList levelDecay_;

        //- Output directory base
        fileName outDirBase_;

        //- Once-per-run header written?
        mutable bool zetaPlanesFileWritten_;

        //- Field-existence check fires once on first execute(). It
        //  cannot run in read()/constructor because some producing
        //  functionObjects (e.g. grad) only register their output
        //  field on their own first execute().
        mutable bool fieldsChecked_;


    // Private member functions

        //- Build the zeta level distribution internally from dz0_,
        //  stretch_, zMaxAir_, zMaxWater_. Walks geometrically outward
        //  from zeta = 0 on both sides, stopping when |zeta| would
        //  exceed the per-side extent. The level at zeta = 0 is
        //  included exactly once.
        void buildZetaLevels();

        //- Write the level list to postProcessing/<name>/zetaLevels.txt
        //  (master only, once per run). One level per line, with the
        //  index, zeta value, assigned Nxy and decay factor D_k.
        //  Convenient for plotting the level distribution and the
        //  per-level resolution.
        void writeZetaLevelsFile() const;

        //- Write zetaPlanes.txt (master only, once per run). Records the
        //  spec the postprocess pipeline needs to read each binary file:
        //  per-level Nxy, zeta value, decay factor, level index, kp,
        //  fields list, and binary layout convention.
        void writeZetaPlanesFile() const;

        //- Resolve which tier a |zeta| value belongs to. Returns the
        //  tier's Nxy, or -1 if the value falls outside every tier.
        //- Lookup tier Nxy for a signed zeta. Picks tiersAir_ if
        //  zeta > 0, tiersWater_ if zeta < 0 (and zeta = 0 uses
        //  whichever list has the smallest innermost tier).
        label tierNxy(scalar zeta) const;

        //- Resample the surfaceElevation grid (Nx_src x Ny_src) onto the
        //  sampler's per-level (nxy x nxy) grid via nearest-neighbor.
        void resampleEtaToLevel
        (
            const List<scalar>& etaSrc,
            label NxSrc, label NySrc,
            scalar dlSrc, scalar xSrc0, scalar ySrc0,
            label nxy,
            List<scalar>& etaOut
        ) const;

        //- Disallow copy / assign
        zetaSampler(const zetaSampler&) = delete;
        void operator=(const zetaSampler&) = delete;


public:

    //- Runtime type information
    TypeName("zetaSampler");


    // Constructors

        zetaSampler
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~zetaSampler() = default;


    // Member Functions

        virtual bool read(const dictionary& dict);
        virtual bool execute();
        virtual bool write();
};


} // End namespace functionObjects
} // End namespace Foam

#endif


