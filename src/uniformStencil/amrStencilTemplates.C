/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2026 Asim Onder
     \\/     M anipulation  |
-------------------------------------------------------------------------------
\*---------------------------------------------------------------------------*/

#include <queue>
#include <cmath>
#include "indexedOctree.H"
#include "treeDataPoint.H"
#include "treeBoundBox.H"
#include "Random.H"
#include "labelIOList.H"

// * * * * * * * * * * * * Private helper * * * * * * * * * * * * * * * * * //

static bool bboxOverlaps
(
    const Foam::point& eMin,
    const Foam::point& eMax,
    const Foam::point& dMin,
    const Foam::point& dMax,
    const Foam::Vector<bool>& isCyclic,
    const Foam::vector& L
)
{
    for (Foam::direction d = 0; d < 3; d++)
    {
        bool dimOK = false;
        if (eMin[d] <= dMax[d] && eMax[d] >= dMin[d])
            dimOK = true;

        if (!dimOK && isCyclic[d] && L[d] > Foam::SMALL)
        {
            if ((eMin[d]+L[d]) <= dMax[d] && (eMax[d]+L[d]) >= dMin[d])
                dimOK = true;
            if ((eMin[d]-L[d]) <= dMax[d] && (eMax[d]-L[d]) >= dMin[d])
                dimOK = true;
        }
        if (!dimOK) return false;
    }
    return true;
}


static bool cellInsideBBox
(
    const Foam::point& cc,
    const Foam::point& bMin,
    const Foam::point& bMax,
    const Foam::Vector<bool>& isCyclic,
    const Foam::vector& L
)
{
    for (Foam::direction d = 0; d < 3; d++)
    {
        bool dimIn = false;
        if (cc[d] >= bMin[d] && cc[d] <= bMax[d])
            dimIn = true;

        if (!dimIn && isCyclic[d] && L[d] > Foam::SMALL)
        {
            if ((cc[d]+L[d]) >= bMin[d] && (cc[d]+L[d]) <= bMax[d])
                dimIn = true;
            if ((cc[d]-L[d]) >= bMin[d] && (cc[d]-L[d]) <= bMax[d])
                dimIn = true;
        }
        if (!dimIn) return false;
    }
    return true;
}


static Foam::Vector<Foam::label> spatialKey
(
    const Foam::point& p,
    Foam::scalar invBin
)
{
    return Foam::Vector<Foam::label>
    (
        Foam::label(std::floor(p.x() * invBin)),
        Foam::label(std::floor(p.y() * invBin)),
        Foam::label(std::floor(p.z() * invBin))
    );
}

static Foam::label packKey(const Foam::Vector<Foam::label>& k)
{
    return k.x() + 99991 * k.y() + 999979 * k.z();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
template<typename Type>
void Foam::amrStencil::getZoneField
(
    const boolList& zone,
    Map<Type>& phiMap,
    const GeometricField<Type, fvPatchField, volMesh>& phi
)
{
    const bool enableTimeLogs = true;

    updateIfNeeded();

    if (!procMapsBuilt_)
    {
        buildProcMaps();
    }

    const label nLayers = max(max(iMax_, jMax_), kMax_);

    Map<DynamicList<label>> spatialHash;
    labelList  remoteGIds;
    pointField remoteCCs;
    Map<Type>  gIdToValue;
    scalar binSize = 0;
    scalar invBin = 0;

    label nBandCells = 0;
    forAll(zone, ci) if (zone[ci]) nBandCells++;

    // ===============================================================
    // Communication phase: dispatch to selected algorithm
    // ===============================================================
    if (Pstream::parRun())
    {
        if (commAlgorithm_ == "allToAll")
        {
            #include "commAllToAll.H"
        }
        else if (commAlgorithm_ == "faceOnly")
        {
            #include "commFaceOnly.H"
        }
        else if (commAlgorithm_ == "faceRelay2Hop")
        {
            #include "commFaceRelay2Hop.H"
        }
        else if (commAlgorithm_ == "bboxDirect")
        {
            #include "commBBoxDirect.H"
        }
        else
        {
            FatalErrorInFunction
                << "Unknown commAlgorithm: " << commAlgorithm_ << nl
                << "Available: allToAll, faceOnly, faceRelay2Hop, "
                << "bboxDirect"
                << abort(FatalError);
        }
    }

    // ===============================================================
    // Common topology build, local collect, spatial-hash fill
    // ===============================================================
    #include "commonTopologyFill.H"

    // ===============================================================
    // Optional missing-stencil rescue phase 
    // ===============================================================
    if (Pstream::parRun() && rescueMissingStencils_)
    {
        #include "commRescueMissingStencils.H"
    }
}
