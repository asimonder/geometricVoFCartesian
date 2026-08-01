# amrStencil: AMR-Compatible Height Function Stencil for OpenFOAM

## Class Hierarchy

```
uniformStencil                 (refactored: private -> protected data)
  |  A_, N_, NxS_, NyS_, NzS_, iMax_, jMax_, kMax_, is2D_, isEmpty_, dx_
  |  a2(), a3()
  |  calcYoungNormal(), calcCCDNormal()
  |  flip(), transpose(), getSmallStencil(), getStencil()
  |  setStencil(phiIJK, ijk3)     <-- IJK-based, uses ijkZone
  |
  +-- amrStencil                  (new: inherits everything above)
        setStencil(phiIJK, seedCelli)  <-- cached topology lookup
        getZoneField(zone, phiMap, phi)  <-- builds + caches topology, parallel exchange
        update()                         <-- call after mesh change
```

## Caching Architecture

Three levels of cached data, all invalidated by `update()`:

| Cache | Built when | Cost | Reused |
|-------|-----------|------|--------|
| `neighbourProcs_` | Constructor / `update()` | Trivial | Always |
| `procPatchRemoteGId_`, `procPatchRemoteCC_` | Constructor / `update()` | One proc exchange | Every `getZoneField` call |
| `stencilTopology_[celli]` | First `getZoneField` call for that cell | One BFS per cell (~63 cells visited) | Every subsequent `setStencil` + `getZoneField` |

**After the first time step on a given mesh**, all topology is cached. `setStencil` becomes a pure array lookup (loop over 63 entries, read from `phiIJK` map). `getZoneField` skips BFS entirely and only scans cached topology to identify remote cells for exchange. Cost is equivalent to the uniform grid version.

## When to call update()

```cpp
if (mesh.changing())   // or mesh.topoChanging() for AMR specifically
{
    stencil.update();
}
```

This clears all caches and rebuilds processor boundary maps. The per-cell topology cache rebuilds lazily on the next `getZoneField` call.

## Usage

```cpp
// Construct once (e.g. in solver setup)
amrStencil stencil(mesh, 3, 3, 7, dx, dy, dz,
                   isEmpty, isCyclic, Lx, Ly, Lz);

// Each time step:
if (mesh.changing())
{
    stencil.update();
}

Map<scalar> alphaMap;
stencil.getZoneField(interfaceCells, alphaMap, alpha);

forAll(interfaceCells, celli)
{
    if (interfaceCells[celli])
    {
        stencil.setStencil(alphaMap, celli);    // cached lookup
        vector n = stencil.calcYoungNormal();   // inherited
        // ... height function code unchanged
    }
}
```

## Changes to uniformStencil

Backward-compatible modifications:

1. `const ijkZone& ijkMesh_` (private) -> `const ijkZone* ijkMeshPtr_` (private, nullable)
2. New protected section: `A_`, `N_`, `NxS_/NyS_/NzS_`, `iMax_/jMax_/kMax_`, `is2D_`, `isEmpty_`, `dx_`, and `SubclassTag` constructor
3. `calcYoungNormal`/`calcCCDNormal` use member `dx_`/`isEmpty_` instead of `ijkMesh_.dx()`/`ijkMesh_.isEmpty()` (values are identical, set from `ijkMesh` at construction)
4. `virtual` destructor added
5. All public signatures preserved exactly

## Known limitations

1. **Single processor hop**: BFS cannot reach cells more than one processor boundary away. Ensure processor subdomains are wider than `kMax * dz`.

2. **dx must match finest AMR level**: Pass the interfacial band dx explicitly.

3. **Memory**: The topology cache stores one `labelList(63)` per interface cell. For 10,000 interface cells that is ~5 KB. Negligible.

4. **No symmetry mirroring**: The BFS-based `setStencil` does not implement `symXIn`/`symXOut` boundary reflection. Can be added if needed.

## Files

| File | Status | Notes |
|------|--------|-------|
| `uniformStencil.H` | Modified | Protected members, SubclassTag ctor, virtual dtor |
| `uniformStencil.C` | Modified | dx_/isEmpty_ members, ijkMeshPtr_ |
| `amrStencil.H` | New | Add to includes |
| `amrStencil.C` | New | Add `uniformStencil/amrStencil.C` to `Make/files` |
| `amrStencilTemplates.C` | New | Included via `#ifdef NoRepository` |
