/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2026 Asim Onder
     \\/     M anipulation  |
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
    Foam::splitAdvectorLim

Description
    Direction-split geometric VOF advection after Weymouth & Yue
    (JCP 229, 2010), stage 1: volume fraction only, 2D uniform
    Cartesian meshes.

    Per time step the two active coordinate directions are swept
    sequentially, order alternating with step parity. Each sweep
    computes geometric liquid fluxes through the faces normal to the
    sweep direction from the donor cell's PLIC reconstruction (normals
    from the reconstructionSchemes registry, e.g. mycCartesian or
    plicRDF), and applies the Weymouth-Yue compression term

        alpha_i += cH_i * (dt/V) * sum_d-faces(phi_f_outward)

    with cH_i = H(alpha_i^n - 1/2) frozen at the step start. Because
    the projected flux phi is discretely solenoidal, the compression
    contributions cancel exactly over the full sweep cycle, so the net
    step update is in conservative flux form and the accumulated
    liquid face flux alphaPhi() is exactly consistent with the alpha
    update. Boundedness holds for face Courant numbers up to 1/2
    without clipping; excursions beyond round-off are fatal by design,
    not repaired.

    Parallel and cyclic boundaries are handled through the coupled
    fvPatchField mechanism (patchNeighbourField of alpha and of the
    reconstructed normal); no zone or stencil communication is used.

    Scope guards (stage 1): 2D (one empty direction), uniform
    orthogonal Cartesian mesh, static mesh, no sources. Violations are
    fatal at construction or at the first advect() call.

SourceFiles
    splitAdvectorLim.C

\*---------------------------------------------------------------------------*/

#ifndef splitAdvectorLim_H
#define splitAdvectorLim_H

#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "reconstructionSchemes.H"
#include "geometricAdvector.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class splitAdvectorLim Declaration
\*---------------------------------------------------------------------------*/

class splitAdvectorLim
:
    public geometricAdvector
{
    // Private data

        volScalarField& alpha1_;
        const surfaceScalarField& phi_;
        const volVectorField& U_;
        const fvMesh& mesh_;
        const dictionary dict_;

        //- Interface reconstruction (registry-selected, e.g. mycCartesian)
        autoPtr<reconstructionSchemes> surfPtr_;

        //- Accumulated liquid volume flux of the completed step [m3/s]
        surfaceScalarField alphaPhi_;

        //- Accumulated total momentum face flux of the completed step
        //  (rho-weighted, both phases summed) on zone faces
        //  [kg m / s^2]; zero elsewhere. Lets a solver stitch an
        //  outside-zone convection operator to the same faces.
        surfaceVectorField qPhi_;

        //- Active sweep directions (component indices) and their count
        //  (2 in 2D with one empty direction, 3 in 3D)
        FixedList<label, 3> sweepDirs_;
        label nSweeps_;

        //- Zone restriction (ijkDict cellSet contract; "domain" or a
        //  missing ijkDict means every cell). ALPHA sweeps and
        //  updates act on zone cells only (PLIC needs the uniform
        //  lattice); MOMENTUM is advected on every axis-aligned face
        //  of the mesh with the same donor + (1-Co)/2 offset
        //  construction, which is the per-phase Lax-Wendroff flux
        //  and needs no lattice. The interface is confined to the
        //  zone interior by the margin guard.
        boolList inZone_;
        boolList inMargin_;
        bool wholeDomain_;

        //- Dynamic finest-level band mode (ijkDict entry dynamicBand):
        //  the zone is every cell at the minimum cell volume, rebuilt
        //  whenever the cell count changes (AMR); margins then sit at
        //  refinement-level jumps
        bool dynamicBand_;

        //- Cell count at the last zone build (rebuild trigger)
        label nCellsAtZoneBuild_;

        //- Frozen Heaviside marker of the current step
        scalarField cH_;

        //- Frozen per-phase compression values of the current step
        //  (stage 2): tc1 = q1/alpha where alpha > 1/2 else 0;
        //  tc2 = q2/(1-alpha) where alpha < 1/2 else 0
        vectorField tc1_, tc2_;

        //- Wisp clamp tolerance (round-off repairs only; larger
        //  excursions are fatal)
        const scalar clampTol_;
        const scalar fatalTol_;

        //- Counters (step)
        label nClamped_;
        scalar maxExcursion_;


    // Private Member Functions

        //- Liquid area fraction of the slab [xa,xb] x [0,1] of a unit
        //  cell cut by a line with inward L1-normalized normal
        //  (mx, mt), mx >= 0, mt >= 0, liquid on n.x <= ai, where ai
        //  is set so the total liquid fraction equals a.
        static scalar slabFraction
        (
            const scalar mx,
            const scalar mt,
            const scalar a,
            const scalar xa,
            const scalar xb
        );

        //- Liquid fraction of the unit cube below the plane
        //  m1 x + m2 y + m3 z = ai, with mk >= 0, m1+m2+m3 = 1
        //  (Scardovelli-Zaleski forward problem, symmetric wrapper)
        static scalar boxFraction3
        (
            const scalar m1,
            const scalar m2,
            const scalar m3,
            const scalar ai
        );

        //- Inverse problem: intercept ai such that boxFraction3 = a
        //  (bracketed Newton on the monotone forward function)
        static scalar boxIntercept3
        (
            const scalar m1,
            const scalar m2,
            const scalar m3,
            const scalar a
        );

        //- Geometric donor fraction for one face: liquid fraction of
        //  the volume swept through the face during dt. d = sweep
        //  component, sOut = +1 if the flux leaves the donor through
        //  its +d face, -1 through its -d face, w = swept width as a
        //  fraction of the donor cell, aD = donor alpha, nD = donor
        //  interface normal (any magnitude/orientation).
        scalar donorFraction
        (
            const label d,
            const scalar sOut,
            const scalar w,
            const scalar aD,
            const vector& nD
        ) const;

        //- Build zone membership and the confinement margin
        void buildZone();

        //- Zone guards: uniform cell volume and axis-aligned faces
        void checkZone();

        //- Common step driver (null q pointers = alpha only)
        void advanceStep(volVectorField* q1, volVectorField* q2);

        //- One directional sweep: update alpha1 internal field,
        //  accumulate alphaPhi_, return max face Courant seen.
        //  When q1/q2 are non-null, the phase momenta are advected in
        //  the same sweep with the same geometric fluxes (stage 2):
        //  face values from majority-phase concentration gradients
        //  with the (1-Co) half-cell offset, compression from the
        //  frozen tc fields (Weymouth-Yue consistency for tracers).
        scalar sweep
        (
            const label d,
            const scalar dt,
            volVectorField* q1,
            volVectorField* q2
        );

        //- Gather along-d neighbour values of a field into plus/minus
        //  cell arrays (coupled patches via patchNeighbourField;
        //  physical boundaries copy the owner value)
        template<class Type>
        void gatherNeighbours
        (
            const label d,
            const GeometricField<Type, fvPatchField, volMesh>& f,
            Field<Type>& fMinus,
            Field<Type>& fPlus
        ) const;

        //- Majority-phase concentration gradient along d (unlimited
        //  centered where both neighbours carry the phase, one-sided
        //  otherwise, zero where the donor does not carry it)
        static vector concGradient
        (
            const scalar cM, const scalar cC, const scalar cP,
            const vector& tM, const vector& tC, const vector& tP,
            const bool inverse
        );

        //- No copy construct / assignment
        splitAdvectorLim(const splitAdvectorLim&) = delete;
        void operator=(const splitAdvectorLim&) = delete;


public:

    //- Runtime type information
    TypeName("splitAdvectorLim");


    // Constructors

        //- Construct from alpha, flux, velocity; reconstruction and
        //  controls read from fvSolution "alpha.*" (same convention
        //  as isoAdvection)
        splitAdvectorLim
        (
            volScalarField& alpha1,
            const surfaceScalarField& phi,
            const volVectorField& U
        );


    //- Destructor
    virtual ~splitAdvectorLim() = default;


    // Member Functions

        //- Advance alpha1 by one time step
        virtual void advect();

        //- Advance alpha1 AND the phase momenta q1 = alpha rho1 u,
        //  q2 = (1-alpha) rho2 u consistently (stage 2). The caller
        //  recovers u = (q1+q2)/rho(alpha) afterwards.
        virtual void advect(volVectorField& q1, volVectorField& q2);

        //- Accumulated liquid volume face flux of the completed step
        virtual const surfaceScalarField& alphaPhi() const
        {
            return alphaPhi_;
        }

        //- Total momentum face flux of the completed step (zone faces)
        const surfaceVectorField& qPhi() const
        {
            return qPhi_;
        }

        //- Zone membership per cell
        const boolList& inZone() const
        {
            return inZone_;
        }

        //- Whole-domain zone?
        bool wholeDomain() const
        {
            return wholeDomain_;
        }

        //- Access the reconstruction
        virtual reconstructionSchemes& surf()
        {
            return *surfPtr_;
        }

        //- Step diagnostics
        virtual label nClamped() const { return nClamped_; }
        virtual scalar maxExcursion() const { return maxExcursion_; }
};


} // End namespace Foam

#endif

// ************************************************************************* //
