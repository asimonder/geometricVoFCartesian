interWYConsistentCenteredJumpFoam
=================================

Fork of interWYConsistentCenteredFoam (2026-07-25), created on his
instruction as a solver separate from interWYConsistentCenteredGFMFoam,
which another session is developing.  Nothing outside this directory is
touched, and neither is geometricVoFCartesian/src.

Parent: applications/solvers/interWYConsistentCenteredFoam, the binary
that ran wave/N0512_Co0.10_cons and carries the recorded +1.252e-4 total
energy loss over 10 s at Co 0.10.

Discipline: every source file has a .pristine copy.  Changes are applied
by applyJump*.py, which restores from pristine first, so edits never
stack, and which asserts every replacement so a silent no-op cannot pass.


Design basis measured on N0512_Co0.10_cons at t = 6 s
-----------------------------------------------------

1.  The gravity term needs no correction for a single-valued interface.
    Per column the telescoped density jump is -998.8000 exactly, and the
    density-gradient-weighted height equals the interface height from the
    column volume fraction to printed precision, with no dependence on
    the sub-cell interface position.  This is an Abel summation identity.
    Relocating the buoyancy jump therefore corrects nothing here, while
    perturbing the exchange between kinetic and potential energy by 33%
    of that term, against a total drift budget of 1.252e-4 over 22667
    steps.  The solver leaves the gravity term alone.

2.  The face density in the pressure equation is the one quantity with a
    defensible error.  The parent uses the harmonic average of the two
    mixture cell densities, which is the one-sided expression with the
    interface pinned to the face midpoint.  Against the true interface
    position the parent coefficient is wrong by a median factor 1.99 and
    by more than 10 on 43.8% of the crossed faces.

3.  rho_f cancels out of rho_f a_f, so the face density cannot act on the
    kinetic-to-potential exchange.  Its effect is confined to the
    pressure solution and to the cell correction.


Step 1, implemented
-------------------

    1/rho_f = lambda/rho_s(P) + (1 - lambda)/rho_s(N)
    lambda  = psi_P/(psi_P - psi_N)

psi is the reconstructed distance function, built each step from the
advector's PLIC planes through reconstructedDistanceFunction, the same
route his rdfCurvature model uses, so no reconstruction scheme is
required to supply it.  Measured sign convention: psi > 0 in phase 1.
A face is treated when both cells are in the band and psi changes sign.
Coupled patches go through syncTools::swapBoundaryCellList and are
treated exactly like internal faces.

Dictionary entries, PIMPLE dict, all defaulting to the parent path:

    jumpFaceDensity     false   one-sided face density on or off
    jumpReconstruct     true    re-reconstruct after the advection so
                                psi and rho describe the same interface
    jumpReport          true    periodic face count and rho_f range
    jumpReportInterval  100


Gates
-----

Phase 0, fork identical to the parent.  N32 wave case, serial, 2201
numerical log lines and the three monitor files byte identical, md5
905cf52cfb84f03bdf3e57ea4569fd9b on both.  PASSED.

Step 1 switch off, identical to the parent.  Same case, monitors byte
identical, md5 ee3e904d2e43591534d2769c55f2adba on both.  PASSED.

Step 1 switch on, N32, 0 to 0.28 s, relative total energy loss:

    parent                  +1.611e-2
    one-sided, reconstruct  +7.485e-2      factor 4.65
    one-sided, no recon     +7.548e-2      factor 4.69

FAILED at this resolution.  The extra reconstruction is not the cause,
the two variants agree to 0.8%.  107 to 108 faces carry the one-sided
value on the 32 x 32 mesh and rho_f spans the full [1.2, 1000] range,
so the treatment is active as intended.

Well-balancing is unaffected by construction.  At rest the required
pressure jump satisfies snGrad(p_rgh) = -ghf snGrad(rho) independently of
rho_f, because the same rho_f divides the gravity term and the pressure
gradient.  The N32 result is a property of the coefficient, not a broken
balance.

Conjecture for the dissipation.  The one-sided coefficient describes a
sharp sub-cell density distribution, while the momentum recovery
u = (q1 + q2)/rho and the cell correction hf/rho_cell both use the
mixture density of the same cell.  In a mixed cell the projection and the
momentum recovery then describe different objects.  That is the pairing
the parent's harmonic-plus-cell-density combination respects.  Untested.

Open: resolution dependence.  The band is 2 cells of 32 here against 2 of
512 in production, so N32 is the worst case for the coefficient spread.
The measured spread at N512 is still median 1.99 with p90 79.9, so it
does not vanish, but the energy penalty may.


Step 2, advected face density: DIVERGES
---------------------------------------

rho_f = rho1 alpha_f + rho2 (1 - alpha_f), alpha_f = phi_alpha/phi clipped
to [0,1] and blended onto interpolate(alpha1) where phi crosses zero.
Switch advectedFaceDensity, default false.  N32 wave case: max|phi|
reaches 2.8e14 by the first report and the run dies.  The volume fraction
and the velocity residuals stay at their initial values while phi runs
away, so the divergence is in the projection, not in the transport.

Isolation run with advectedPhiTol 1e30, which makes alpha_f exactly
interpolate(alpha1), diverges the same way.  The instability is therefore
the arithmetic alpha-weighted face density itself, not the alpha-flux
ratio or its regularisation.

Both candidate coefficients now fail.  The one-sided value is stable but
4.65x more dissipative at N32.  The arithmetic value diverges.  The
parent's harmonic average of the mixture cell densities lies between them
and is the only stable, low-dissipation choice found.  The face density in
the Poisson coefficient is not where further effort should go.


Step 3, cell correction from the face acceleration: INJECTS
-----------------------------------------------------------

gC = avg(gf) instead of avg(rho_f gf)/rho_cell.  Switch
cellCorrectionAcceleration, default false, switch-off regression byte
identical to the parent on all three monitors.

N32, 0 to 0.28 s, relative total energy:

    parent            +1.611e-2   loss
    gC acceleration   -5.828e-1   GAIN
    maxUwater at T:   0.268 parent, 0.464 variant

Mechanism.  At rest gf = 0 on every face, so this is not a balance
failure.  In motion the interface faces carry a residual gf, the part of
the face force the exactly projected flux does not cancel, and it is
amplified by 1/rho_f with rho_f harmonic and of order 2.4 kg/m3.  The
parent multiplies that residual by rho_f and divides by rho_cell, which
attenuates it by rho_f/rho_cell, of order 5e-3.  Averaging the
acceleration removes that attenuation and the residual drives the flow.

Consequence.  The parent's combination of harmonic face density, force
averaging and mixture cell density is a tightly coupled set.  All three
perturbations tried in this solver degrade it: the one-sided coefficient
by 4.65x in dissipation, the arithmetic coefficient into divergence, the
acceleration average into 58% energy gain.  The pressure-velocity stage is
closed to interface-sharpening changes.

It also locates the interface noise.  It is the residual gf at interface
faces, the mismatch between the exactly projected face flux and the
collocated cell velocity, attenuated by 5e-3 and leaking through.  It is a
property of the approximate projection, not of any missing jump condition.
