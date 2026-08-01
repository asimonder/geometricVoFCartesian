# surfaceElevation Function Object

Online extraction of $\eta(x, y, t)$ via local height-function column integration on the AMR mesh, implementing the algorithm described in [[dns_doc#4.1.1 Surface elevation extraction]].

## Outline

- [[#1. What it does]]
- [[#2. How the column stencil is realised]]
- [[#3. Failure modes and counters]]
- [[#4. Output layout]]
- [[#5. Statistics output]]
- [[#5a. Per-time validation diagnostics]]
- [[#6. Build]]
- [[#7. Wiring into a case]]
- [[#8. Limitations]]
- [[#9. Open work]]

## 1. What it does

For every cell with $\mathrm{interfaceTol} < \alpha < 1 - \mathrm{interfaceTol}$, the function object builds a $1 \times 1 \times 7$ vertical column of cells centred on the interface cell, verifies that the bottom cell is fully wet and the top cell is fully dry, and computes

$$
\eta_{\text{local}} \,=\, \Delta z \sum_{k=-3}^{+3} \alpha_k
\,, \qquad
\eta_{\text{abs}} \,=\, z_{\text{seed}} - 3.5\,\Delta z \,+\, \eta_{\text{local}} \,,
$$

with $z_{\text{seed}}$ the cell-centre $z$ of the interface cell. The wave excursion is $\eta_{\text{abs}} - z_0$ where $z_0$ is the still-water reference. The construction is exact on a flat interface and second-order accurate on a curved one, in line with the doc.

## 2. How the column stencil is realised

The 1×1×7 column is the existing `amrStencil` constructed with `iMax = 0`, `jMax = 0`, `kMax = 3`. The BFS in `amrStencil::buildStencilTopology` rejects any neighbour outside the column via `withinStencil()`, so the BFS walks straight up and down without ever exploring lateral cells. The cross-processor traversal that already works for height-function curvature applies here unchanged. No new amrStencil functionality, no new communication code.

The construction mirrors `heightFunctionAMR` for the cyclic / empty / boundary detection so behaviour at boundaries is identical to the curvature stencil. Cells within `kMax` layers of a wall, inlet, outlet or symmetry boundary are skipped (they cannot have a complete column). Cyclic and processor boundaries are not flagged because the stencil wraps through them.

## 3. Failure modes and counters

Five categories are counted at every output and printed both to the console and to the per-time `eta_summary.dat`:

| counter | meaning |
|---|---|
| `nIncomplete` | BFS did not reach all 7 cells (typically: column crossed a processor boundary in $z$ within reach of the seed) |
| `nMixedRefinement` | one or more column cells differ from the seed in volume by more than 2x (column straddles AMR refinement transition) |
| `nBadBottom` | $\alpha$ of the bottom cell of the column is less than $1 - \mathrm{wetTol}$, i.e. the column does not extend deep enough to reach water |
| `nBadTop` | $\alpha$ of the top cell of the column exceeds $\mathrm{dryTol}$, i.e. does not reach air |
| `nMultipleInterf` | subset of `nBadBottom + nBadTop` where the failing end-cell is itself an interface cell, signalling a microbreaking / overhang event in that column |

Each failed column is silently dropped from the output rather than written with a sentinel, on the grounds that downstream Python postprocess prefers missing entries to flagged ones.

## 4. Output layout

Each rank writes one ASCII file per output time at

```
<case>/postProcessing/surfaceElevation/<time>/eta_proc<N>.dat
```

with one row per valid column:

```
x  y  eta_abs  eta_excursion
```

Master writes `eta_summary.dat` in the same directory containing the global counters and the validation diagnostics.

For a regular-grid representation (as required by §4.2 spectra), reading the rank files and binning on $(x, y)$ to the AMR Lvl1 cell pitch is trivial because the interfacial band is uniform: `i = round((x - x_min) / dx)`, `j = round((y - y_min) / dy)`. A 12-line numpy snippet does the merge.

## 5. Statistics output

Master writes a single growing time-series file at

```
<case>/postProcessing/surfaceElevation/<startTime>/surfaceElevation.dat
```

via the `writeFile` base, with one row per output time. Columns:

```
# time   sigmaEta   Hs   etaMin   etaMax   HmaxOverHs   skew   kurt
```

`sigmaEta` is the standard deviation of the wave excursion. `Hs = 4·sigmaEta` is the significant wave height in the engineering convention. `HmaxOverHs = (etaMax - etaMin) / Hs` is the rogue-wave criterion; values above 2 indicate exceptional events. `skew` and `kurt` are the third and fourth standardised moments and test the Gaussianity assumption underlying linear wave theory; `kurt` is raw (not excess), subtract 3 if excess is wanted.

This column set was settled informally and has not been cross-checked against the wind-wave DNS literature. Wu et al. 2022, Li & Shen 2022 and Laxague et al. 2018 emphasise mean square slope $\langle |\nabla\eta|^2 \rangle$ and the streamwise/spanwise slope-variance components $\langle\eta_x^2\rangle$, $\langle\eta_y^2\rangle$ as the primary growth and anisotropy diagnostics; these are not yet computed. Adding them requires either a second `amrStencil` with `iMax=jMax=1, kMax=0` for horizontal neighbours, or a refactor that extracts a shared height field from `heightFunctionAMR` for both consumers (curvature and elevation).

## 5a. Per-time validation diagnostics

The per-time `eta_summary.dat` reports the failure-counter breakdown (`nIncomplete`, `nMixedRefinement`, `nBadBottom`, `nBadTop`, `nBoundarySkip`, `nMultipleInterf`). It is the file to inspect when validating a run.

## 6. Build

The function object is part of the existing `libgeometricVoFCartesian` library, alongside `heightDistanceFunction` and the rest of the package. Source goes under `src/functionObjects/surfaceElevation/`:

```
geometricVoFCartesian/
└── src/
    ├── ijkZone/
    ├── uniformStencil/
    ├── interfaceForces/
    ├── functionObjects/
    │   ├── heightDistanceFunction/
    │   └── surfaceElevation/                <- new
    │       ├── surfaceElevation.{H,C}
    │       └── README.md
    └── Make/
        ├── files
        └── options
```

There is no separate `Make/` subdirectory under the new folder. Add one line to `src/Make/files` next to the existing `heightDistanceFunction` entry:

```
functionObjects/heightDistanceFunction/heightDistanceFunction.C
functionObjects/surfaceElevation/surfaceElevation.C
```

`src/Make/options` does not need editing; the existing includes already cover `finiteVolume`, `fieldFunctionObjects`, `dynamicMesh`, etc., and `amrStencil.H` / `heightFunctionAMR.H` are inside the same library and visible without further `-I`.

Rebuild from `src/`:

```bash
cd $WM_PROJECT_USER_DIR/geometricVoFCartesian/src
wmake libso
```

The function object is then inside `$FOAM_USER_LIBBIN/libgeometricVoFCartesian.so`.

## 7. Wiring into a case

In `system/controlDict`:

```
functions
{
    surfaceElevation
    {
        type            surfaceElevation;
        libs            ("libgeometricVoFCartesian.so");

        executeControl  adjustableRunTime;
        executeInterval 0.05;       // T_p / 8 for Case 1
        writeControl    adjustableRunTime;
        writeInterval   0.05;       // same cadence as executeInterval

        alpha           alpha.water;
        z0              0.0;        // still-water reference

        kMax            3;          // 7-cell column
        wetTol          1.0e-6;
        dryTol          1.0e-6;
        interfaceTol    1.0e-12;
    }
}
```

The function object reads alpha by name from the registry, so it works in both the live solver run and the foamRun-style replay path.

## 8. Limitations

The BFS in `amrStencil` does not walk multi-hop through processor boundaries. Columns whose seed cell lies within $\mathrm{kMax}$ cells of a $z$-direction processor decomposition boundary will be incomplete and dropped (counted in `nIncomplete`). This is the same constraint as `heightFunctionAMR` and was acceptable there for the same reason: lateral $(x, y)$ decomposition keeps full vertical columns inside one rank, and column failures only occur when $z$-decomposition cuts the AMR-refined band. For the production decomposition this is rare or zero.

Microbreaking columns where two interface cells stack within the 7-cell window are dropped via `nBadBottom` or `nBadTop` and counted separately in `nMultipleInterf`. A widened column (`kMax = 4` or `5`) recovers the deeper of two stacked interfaces if the band is thick enough; the trade-off is more `nMixedRefinement` failures at the band edge, where AMR transitions to coarser cells. The default `kMax = 3` matches the doc's specification and the typical 8-10-cell vertical AMR band.

## 9. Open work

This file is shipped as a working baseline. It does the column extraction, the integrity checks, and the zeroth- through fourth-order moments of $\eta$. It does not yet do anything that requires horizontal neighbour information.

Items left for the next iteration, in roughly the order they should be tackled:

The first step is the static-water sanity test. Set up a case with uniform $\alpha=1$ below $z=0$ and $\alpha=0$ above, no flow, run for one timestep with the function object enabled. Every valid column should report $\eta_{\mathrm{abs}} = 0$ to machine precision. `sigmaEta`, `etaMin`, `etaMax`, `skew`, `kurt` should all be at machine zero. If they are not, the column algorithm or the boundary handling is broken and no further work is meaningful until that is fixed.

The slope statistics ($\langle\eta_x^2\rangle$, $\langle\eta_y^2\rangle$, $\mathrm{mss} = \langle |\nabla\eta|^2 \rangle$) are the canonical wind-wave DNS diagnostics from the literature and are not yet computed. Two implementation paths were discussed and not pursued: (a) add a second `amrStencil` with `iMax=jMax=1, kMax=0` inside surfaceElevation to fetch lateral neighbours through the existing cross-rank machinery; (b) refactor the column extraction out of `heightFunctionAMR` into a shared height-field class that both curvature and elevation consume, then read the heights and central-difference them. (b) eliminates the redundant column work shared between curvature and elevation but is a larger refactor.

The column set in the time-series file ( `sigmaEta Hs etaMin etaMax HmaxOverHs skew kurt`) was not finalised against the wind-wave DNS literature. `Hs` and `HmaxOverHs` are engineering / extreme-statistics descriptors more typical of ocean engineering than wind-wave growth physics. The wind-wave DNS literature emphasises slope variance and mean square slope, which are not present. The column set should be reviewed once the slope work is done.

The per-rank ASCII files contain the raw $(x, y, \eta)$ data that the §4.2 spectral pass and any other downstream analysis will consume. Whether the per-rank ASCII format is the right output for that pipeline (as opposed to a binned 2D array per rank, or master-only assembled output) is an open design question for whoever builds the spectral diagnostic.
