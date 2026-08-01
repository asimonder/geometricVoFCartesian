/*---------------------------------------------------------------------------*\
            Copyright (c) 2026, Asim Onder
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::functionObjects::surfaceElevation

Group
    grpFieldFunctionObjects

Description
    Online extraction of the air-water interface elevation eta(x, y, t)
    via local height-function column integration on the AMR mesh.

    For each interface cell (a cell with 0 < alpha < 1) the functionObject
    builds a 1 x 1 x 7 vertical column of cells centred on the interface
    cell, verifies that the bottom cell is fully wet and the top cell is
    fully dry, and computes

        eta_local = dz * sum_{k=-3..3} alpha_k

    The absolute interface position is then

        eta_abs = z_seed - 3.5 * dz + eta_local

    with z_seed the cell-centre z of the interface cell. The wave excursion
    is eta_abs - z0 where z0 is the still-water reference.

    The 1 x 1 x 7 stencil is realised by constructing amrStencil with
    iMax = jMax = 0 and kMax = 3. The BFS in amrStencil::buildStencilTopology
    then automatically rejects any neighbour outside the column via
    withinStencil(), and cross-processor vertical traversal works through
    the existing processor / processorCyclic boundary handling.

    Failure modes are caught and counted as diagnostics:
      - incomplete stencil (BFS did not reach all 7 cells)
      - mixed refinement (volume ratio outside [0.5, 2.0])
      - bottom cell not fully wet
      - top cell not fully dry

    Each MPI rank writes its valid columns to
        postProcessing/<name>/<time>/eta_proc<N>.dat
    as ASCII rows of (x  y  eta_abs  eta_excursion). Rank 0 also writes
        postProcessing/<name>/<time>/eta_summary.dat
    with global counts and the volume-conservation and reconstructed
    potential-energy diagnostics described in DNS doc Section 4.1.1.

Usage
    Add to system/controlDict:

    functions
    {
        surfaceElevation
        {
            type            surfaceElevation;
            libs            ("libgeometricVoFCartesian.so");
            executeControl  adjustableRunTime;
            executeInterval 0.05;     // T_p / 8 typical
            writeControl    adjustableRunTime;
            writeInterval   0.05;     // same cadence as executeInterval
            alpha           alpha.water;
            z0              0.0;      // still-water reference
            kMax            3;        // half-stencil in z (default 3)
            wetTol          1e-6;     // (1 - alpha) tol for fully-wet bottom
            dryTol          1e-6;     // alpha tol for fully-dry top
            interfaceTol    1e-12;    // 0 < alpha < 1 detection tol
        }
    }

SourceFiles
    surfaceElevation.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_surfaceElevation_H
#define functionObjects_surfaceElevation_H

#include "fvMeshFunctionObject.H"
#include "volFieldsFwd.H"
#include "Vector.H"
#include "boolList.H"
#include "DynamicList.H"
#include "autoPtr.H"
#include "amrStencil.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                      Class surfaceElevation Declaration
\*---------------------------------------------------------------------------*/

class surfaceElevation
:
    public fvMeshFunctionObject
{
private:

    // Configuration

        //- Name of the volume-fraction field
        word alphaName_;

        //- Still-water reference elevation
        scalar z0_;

        //- Half-stencil height in z (default 3 -> 7-cell column)
        label kMax_;

        //- alpha tolerance for "fully wet": require alpha_bottom > 1 - wetTol
        scalar wetTol_;

        //- alpha tolerance for "fully dry": require alpha_top < dryTol
        scalar dryTol_;

        //- alpha tolerance for "interface cell": interfaceTol < alpha < 1 - interfaceTol
        scalar interfaceTol_;

        //- Output filename for the wave-statistics time series.
        //  Written at the case root, like hamiltonian.dat. Default:
        //  surfaceElevation.dat
        word statsFile_;

        //- Whether to dump the binary 2D eta(x,y) snapshot at every write.
        //  Default: true.
        bool writeEtaField_;


    // Cached grid spec (built once on first write, master-only)

        mutable label  gridNx_;
        mutable label  gridNy_;
        mutable scalar gridDl_;
        mutable scalar gridXmin_;
        mutable scalar gridYmin_;
        mutable bool   gridSpecBuilt_;

        //- Binned eta excursion (eta - z0) on the regular Nx*Ny grid.
        //  Built every execute() on master, broadcast to all ranks, and
        //  stored as a registered IOList<scalar> on the mesh registry
        //  under the name "<this fO name>:etaGrid". Downstream
        //  functionObjects (e.g. zetaSampler) read it via a standard
        //  mesh_.lookupObject<IOList<scalar>>(...) call. NO_WRITE so
        //  it never goes to disk through the standard mechanism.
        //  Indexed [iy * Nx + ix]. Missing cells contain NaN.
        autoPtr<IOList<scalar>> etaGridIOPtr_;

        //- Grid spec stored as a registered IOdictionary on the mesh
        //  registry under "<this fO name>:etaGridSpec", carrying Nx,
        //  Ny, dl, xmin, ymin. NO_WRITE.
        autoPtr<IOdictionary> etaGridSpecPtr_;

        //- Wave-following coordinate zeta = z - eta(x, y, t) defined on
        //  every cell. Updated every execute(). NO_WRITE: not dumped by
        //  surfaceElevation itself; the user can sample it via standard
        //  surfaces / isoSurface functionObjects in controlDict for
        //  visual sanity checks of fields near the interface.
        autoPtr<volScalarField> zetaPtr_;

        //- Finest interface-cell volume (set in execute(), read in write())
        scalar fineV_;


    // Mesh-related state

        //- Empty-direction flags (detected from emptyPolyPatch / wedgePolyPatch)
        Vector<bool> isEmpty_;

        //- Cyclic-direction flags (detected from cyclic / processorCyclic)
        Vector<bool> isCyclic_;

        //- Domain extents
        scalar Lx_;
        scalar Ly_;
        scalar Lz_;

        //- Cell thickness in the empty direction (for 2D)
        scalar emptyDirThickness_;

        //- 2D flag
        bool is2D_;

        //- Cells within kMax of a non-cyclic, non-processor wall boundary
        boolList boundaryCells_;


    // Stencil

        //- Column-only AMR stencil (iMax = jMax = 0, kMax = kMax_)
        autoPtr<amrStencil> stencilPtr_;


    // Extracted data (populated by execute(), flushed by write())

        //- Per-column points: (x, y, eta_abs)
        DynamicList<vector> validPoints_;

        //- Global diagnostic counters (set on master after reductions)
        label gValid_;
        label gIncomplete_;
        label gCoarseLevel_;
        label gBadBottom_;
        label gBadTop_;
        label gBoundarySkip_;
        label gMultipleInterface_;

        //- Wave statistics (area-weighted moments of eta_excursion)
        scalar etaMin_;
        scalar etaMax_;
        scalar etaMean_;
        scalar etaVar_;
        scalar etaStd_;
        scalar etaSkew_;
        scalar etaKurt_;

        //- Slope statistics
        scalar etaXvar_;
        scalar etaYvar_;
        scalar mss_;
        label gValidSlope_;
        label gIncompleteSlope_;


    // Private members

        //- Mark cells within nLayers of a wall / inlet / outlet / symmetry.
        //  These get skipped because the stencil would extend out of the
        //  domain. Cyclic / processor / empty patches are not flagged.
        void markBoundaryCells(label nLayers);

        //- Detect empty / cyclic boundary configuration (mirrors
        //  heightFunctionAMR constructor logic).
        void detectBoundaryConfig();

        //- Write column header to the time-series .dat file. Called once
        //  from write() the first time it opens the file.
        void writeStatsHeader(std::ofstream& os) const;

        //- Write the regular-grid spec to <case>/postProcessing/
        //  surfaceElevation/grid.txt. Called once on master.
        void writeGridSpec() const;

        //- Build the binned eta(x,y) grid on master, broadcast to all
        //  ranks, store in etaGrid_. Called every execute(). Downstream
        //  functionObjects (e.g. zetaSampler) consume the grid via the
        //  public accessors below.
        void buildEtaGrid();

        //- Write etaGrid_ to a binary file at
        //  <case>/postProcessing/surfaceElevation/<time>/eta.bin.
        //  Master-only, called at write cadence.
        void writeEtaBinary() const;

        //- Disallow copy / assign
        surfaceElevation(const surfaceElevation&) = delete;
        void operator=(const surfaceElevation&) = delete;


public:

    //- Public accessors. The eta grid and its spec are exposed via the
    //  mesh objectRegistry under fixed names; downstream functionObjects
    //  (e.g. zetaSampler) read them directly with mesh.lookupObject<>()
    //  rather than calling these accessors. The names below are part of
    //  the public API.

        //- Registry name for the binned eta grid (IOList<scalar>).
        //  Returns "<this fO name>:etaGrid".
        word etaGridName() const { return name() + ":etaGrid"; }

        //- Registry name for the eta grid spec (IOdictionary).
        //  Returns "<this fO name>:etaGridSpec".
        word etaGridSpecName() const { return name() + ":etaGridSpec"; }

        //- Registry name for the zeta volScalarField.
        //  Returns "zeta".
        word zetaName() const { return "zeta"; }

        //- Whether the grid has been built and registered yet.
        bool gridReady() const { return gridSpecBuilt_; }


public:

    //- Runtime type information
    TypeName("surfaceElevation");


    // Constructors

        //- Construct from Time and dictionary
        surfaceElevation
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~surfaceElevation() = default;


    // Member Functions

        //- Read settings (re-callable when controlDict reloads)
        virtual bool read(const dictionary& dict);

        //- Run the column extraction and write per-rank output
        virtual bool execute();

        //- No-op (write is done by execute()).
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
