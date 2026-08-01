interWYConsistentCenteredGFMFoam
Fork of interWYConsistentCenteredFoam (2026-07-25). Ghost Fluid Method
pressure and viscosity treatment at the interface, WY consistent
transport and centered single-projection cycle retained. All GFM
machinery behind a runtime switch, default off = production path bit
for bit. Regenerate from the pristine copy plus the applyGFM*.py
scripts, never stacked hand edits.

Phase 0 (this state): fork, identical to production. Gate: static
  digit identity.
Phase 1: pressure + gravity jump conditions.
  Interface faces: (alpha_P - 0.5)(alpha_N - 0.5) < 0.
  Subface position lambda from the plicRDF reconstructed distance
  function (registered field RDF), not from alpha ratios.
  One-sided beta = 1/rho interpolation + jump source
  H = (rho+ - rho-) g.x_Gamma in the p_rgh Laplacian and in the
  face-force/cell-correction chain (Vukcevic et al. 2017, Eqs. 39-49).
  AMENDMENT (2026-07-25, gate-B finding): at marked faces the cell
  correction delivers the continuous face acceleration beta grad p_d
  to both adjacent cells directly (Eq. 44 discretised as the
  interface-corrected cell gradient, Eqs. 62-65). A force/rho_cell
  delivery there is linearly unstable with gain rho_eff/rho_cell
  (~100/step measured at a wisp row on the N32 flat gate).
  STATUS: built (applyGFMPhase1.py + gfmPressureJump.H); gate A
  byte-identical monitors with the switch off; gate B stable,
  water |Uz|_max 1.6e-20 vs 1.6e-12 production path, same case.
  Gates: hydrostatic rest at machine zero, translating droplet,
  N128 interface noise metric (target: 10x reduction), wave ladder.
Phase 2: viscous jump. COMMITTED SCOPE, not optional.
  Subcell-weighted effective face viscosity (Hwang 2025, Eq. 2.21)
  in the implicit viscous matrix at interface faces, plus the
  explicit jump source J from the full tangential/normal
  decomposition (Hwang Eqs. 2.10-2.13, 2.19-2.20; Sussman 2007
  policy: retain all contributions). Interface geometry (n, tau1,
  tau2) from PLIC / height-function machinery.
  Deferred with evidence: only the viscous term inside the pressure
  jump itself (negligible per Peltonen et al. 2020).
  Gates: two-phase shear-layer suite, wave regression.
Phase 3: surface-tension jump gamma*kappa in [p_d] via the
  height-function curvature at the interface location (Hwang
  Secs. 2.4, 2.7). Gates: static drop, stationary bubble.
