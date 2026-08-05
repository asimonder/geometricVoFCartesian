/*---------------------------------------------------------------------------*\
    Part of interWYGeoConsistentFoam.

    Geometric face density (Liu, Tolle, Bothe, Maric, JCP 493 (2023)
    112426, Eqs. (37)-(39), transplanted from signed distances to the
    PLIC reconstruction): every internal face touching an interface
    cell is cut by that cell's reconstructed plane, and the wetted
    area fraction replaces linear interpolation,

        alpha_f = A_wet/A_f,    rho_f = alpha_f rho1 + (1-alpha_f) rho2.

    The below/above orientation of the cut is calibrated per face from
    the alpha ordering of the two adjacent cells, so no sign convention
    of the reconstruction is assumed. Faces not touching the interface,
    boundary faces, and processor faces keep linear interpolation.
\*---------------------------------------------------------------------------*/

surfaceScalarField alphafGeo(fvc::interpolate(alpha1));
{
    const reconstructionSchemes& surfGeo =
        mesh.lookupObject<reconstructionSchemes>("reconstructionScheme");

    const boolList& icell = surfGeo.interfaceCell();
    const volVectorField& nrmGeo = surfGeo.normal();
    const volVectorField& ctrGeo = surfGeo.centre();

    cutFacePLIC cutter(mesh);

    const labelUList& ownG = mesh.owner();
    const labelUList& neiG = mesh.neighbour();
    const vectorField& SfG = mesh.Sf().primitiveField();
    const scalarField& magSfG = mesh.magSf().primitiveField();
    const vectorField& CG = mesh.C().primitiveField();

    scalarField& afG = alphafGeo.primitiveFieldRef();

    forAll(afG, facei)
    {
        const label o = ownG[facei];
        const label n = neiG[facei];

        scalar sum = 0.0;
        label cnt = 0;

        for (label pass = 0; pass < 2; ++pass)
        {
            const label c = (pass == 0) ? o : n;
            if (!icell[c] || mag(nrmGeo[c]) < VSMALL) continue;

            const vector nh = nrmGeo[c]/mag(nrmGeo[c]);
            const label st = cutter.calcSubFace(facei, nh, ctrGeo[c]);

            scalar fcut;
            if (st == 0)       fcut = mag(cutter.subFaceArea())/magSfG[facei];
            else if (st == -1) fcut = 1.0;   // fully below the plane
            else               fcut = 0.0;   // fully above

            // orientation: the liquid-richer adjacent cell decides
            // whether "below" is the liquid side
            const label rich = (alpha1[o] >= alpha1[n]) ? o : n;
            const bool richBelow = (((CG[rich] - ctrGeo[c]) & nh) < 0);
            sum += richBelow ? fcut : (1.0 - fcut);
            ++cnt;
        }

        if (cnt)
        {
            afG[facei] = min(max(sum/cnt, 0.0), 1.0);
        }
    }

    // instrumentation: face family census and mobility contrast
    label nGeo = 0, nWet = 0;
    scalar afMin = GREAT, afMax = -GREAT;
    forAll(afG, facei)
    {
        const label o = ownG[facei];
        const label n = neiG[facei];
        if (icell[o] || icell[n])
        {
            ++nGeo;
            afMin = min(afMin, afG[facei]);
            afMax = max(afMax, afG[facei]);
            if (afG[facei] > 1e-6) ++nWet;
        }
    }
    Info<< "geoRhof: interface-adjacent faces " << nGeo
        << " wetted " << nWet
        << " alphaf [" << afMin << ", " << afMax << "]" << endl;
}
const surfaceScalarField rhof(alphafGeo*rho1 + (1.0 - alphafGeo)*rho2);
