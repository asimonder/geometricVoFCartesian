#include "cyclicPolyPatch.H"

// interpolateCutCellsToFaces3D.H
// Assumes the following are in scope: alpha1_, K_, Kf_, interfaceCells, interfaceTol_, dl_
// If not already available in the including .C file, add: #include "cyclicPolyPatch.H"

const fvMesh& mesh = alpha1_.mesh();

// 1) Clean up cell-centered curvature where a cell is flagged interface but alpha is effectively pure
forAll(interfaceCells, iCell)
{
    if (interfaceCells[iCell])
    {
        const scalar a = alpha1_[iCell];
        if (a < interfaceTol_ || a > 1.0 - interfaceTol_)
        {
            K_[iCell] = 0.0;
        }
    }
}

// 2) Face interpolation baseline
K_.correctBoundaryConditions();
Kf_ = fvc::interpolate(K_);

// 3) Interior faces: if exactly one side is a cut cell, use the cut-cell curvature
const label nIntFaces = mesh.nInternalFaces();
const labelList& own = mesh.faceOwner();
const labelList& nei = mesh.faceNeighbour();

for (label f = 0; f < nIntFaces; ++f)
{
    const label cO = own[f];
    const label cN = nei[f];

    const scalar aO = alpha1_[cO];
    const scalar aN = alpha1_[cN];

    const bool cutO = (aO > interfaceTol_ && aO < 1.0 - interfaceTol_);
    const bool cutN = (aN > interfaceTol_ && aN < 1.0 - interfaceTol_);

    if (cutO ^ cutN)   // exactly one cut
    {
        Kf_[f] = cutO ? K_[cO] : K_[cN];
    }
    // else: both cut -> keep linear Kf_[f]; both pure -> Kf_[f] should already be ~0
}

// 4) Boundary faces: adjust ONLY true cyclic patches (never processor/wall/etc.)
// If coupled interpolation produced ~0.5*K_owner (neighbor effectively pure), restore K_owner
const scalar invL   = 1.0 / dl_;   // curvature scale ~ 1/Δ
const scalar absTol = 1e-10 * invL; // absolute tol with units of curvature
const scalar relTol = 1e-8;        // relative tol (unitless)

forAll(mesh.boundaryMesh(), iPatch)
{
    const polyPatch& p = mesh.boundaryMesh()[iPatch];

    //if (!isA<cyclicPolyPatch>(p)) continue;  // skip processor, wall, empty, etc.
    if (!mesh.boundaryMesh()[iPatch].coupled()) continue;
    
    auto& Kfb = Kf_.boundaryFieldRef()[iPatch];
    const labelUList& faceCells = p.faceCells();

    forAll(faceCells, iFace)
    {
        const label c    = faceCells[iFace];
        const scalar Kc  = K_[c];
        const scalar Kff = Kfb[iFace];

        if (mag(Kc) > SMALL)
        {
            const scalar thr = max(relTol * mag(Kc), absTol);
            // Detect the "halved owner" case robustly and restore it
            if (mag(Kff - 0.5 * Kc) < thr)
            {
                Kfb[iFace] = Kc;
            }
        }
    }
}
