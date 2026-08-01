/*---------------------------------------------------------------------------*\
    Part of interWYCenteredFoam (see the main file for the full
    description and license). Extracted for the conventional
    fragmented solver layout.
\*---------------------------------------------------------------------------*/

        // momentum is advected on every cell by the sweeps (geometric
        // on the lattice, per-phase Lax-Wendroff in the pure-phase
        // bulk), so the one-fluid velocity is recovered everywhere
        U.primitiveFieldRef() =
            (q1.primitiveField() + q2.primitiveField())
           /rho.primitiveField();
        U.correctBoundaryConditions();
