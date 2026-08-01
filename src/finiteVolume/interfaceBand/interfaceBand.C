/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2026 Asim Onder
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

\*---------------------------------------------------------------------------*/

#include "interfaceBand.H"
#include "FaceCellWave.H"
#include "topoDistanceData.H"
#include "bitSet.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(interfaceBand, 0);
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::word Foam::interfaceBand::findAlphaName(const fvMesh& mesh)
{
    // The two-phase mixture names the field alpha.<firstPhase>, with the
    // phase names declared in transportProperties.
    if (mesh.foundObject<IOdictionary>("transportProperties"))
    {
        const auto& transportDict =
            mesh.lookupObject<IOdictionary>("transportProperties");

        wordList phases;
        if (transportDict.readIfPresent("phases", phases) && phases.size())
        {
            const word name("alpha." + phases[0]);
            if (mesh.foundObject<volScalarField>(name))
            {
                return name;
            }
        }
    }

    // Last resort: common names
    for
    (
        const word& name
      : {"alpha.water", "alpha.phase1", "alpha.liquid", "alpha1"}
    )
    {
        if (mesh.foundObject<volScalarField>(name))
        {
            return name;
        }
    }

    FatalErrorInFunction
        << "Cannot find a volume-fraction field "
        << "(tried alpha.water, alpha.phase1, alpha.liquid, alpha1)." << nl
        << "Set alpha in the interfaceBand sub-dictionary of fvSolution "
        << "to name it explicitly."
        << abort(FatalError);

    return word::null;
}


void Foam::interfaceBand::readControls()
{
    const fvMesh& mesh = band_.mesh();
    const dictionary& solDict = mesh.solutionDict();

    // Legacy flat keys as fallback
    label nLayers =
        solDict.getOrDefault<label>("nBufferLayers", 2);
    scalar threshold =
        solDict.getOrDefault<scalar>("interfaceThreshold", 1e-6);
    word alphaName =
        solDict.getOrDefault<word>("interfaceBandAlpha", word::null);

    // The interfaceBand sub-dictionary wins
    label nRampLayers = 0;
    if (const dictionary* bandDict = solDict.findDict("interfaceBand"))
    {
        nLayers = bandDict->getOrDefault<label>("nLayers", nLayers);
        nRampLayers = bandDict->getOrDefault<label>("nRampLayers", 0);
        threshold = bandDict->getOrDefault<scalar>("threshold", threshold);
        alphaName = bandDict->getOrDefault<word>("alpha", alphaName);
    }

    if (alphaName.empty())
    {
        alphaName = findAlphaName(mesh);
    }

    if
    (
        nLayers != nLayers_
     || nRampLayers != nRampLayers_
     || threshold != threshold_
     || alphaName != alphaName_
    )
    {
        Info<< "interfaceBand: alpha = " << alphaName
            << ", nLayers = " << nLayers
            << ", nRampLayers = " << nRampLayers
            << ", threshold = " << threshold << endl;

        if (nLayers < 1)
        {
            WarningInFunction
                << "nLayers = " << nLayers
                << ": band reduces to mixed cells only; mass-momentum "
                << "consistency margin at the band edge is lost." << endl;
        }
    }

    nLayers_ = nLayers;
    nRampLayers_ = nRampLayers;
    threshold_ = threshold;
    alphaName_ = alphaName;
}


void Foam::interfaceBand::rebuild()
{
    readControls();

    const fvMesh& mesh = band_.mesh();

    const volScalarField& alpha =
        mesh.lookupObject<volScalarField>(alphaName_);

    bitSet isMixed(mesh.nCells());
    forAll(alpha, celli)
    {
        if
        (
            alpha[celli] > threshold_
         && alpha[celli] < (1.0 - threshold_)
        )
        {
            isMixed.set(celli);
        }
    }

    List<topoDistanceData<label>> cellDist(mesh.nCells());

    // topoDistanceData (v2506): updateCell copies, updateFace increments.
    // Seed faces of mixed cells carry distance 1, so a cell k face-hops
    // outside the mixed set receives exactly k; unreached cells stay -1.
    // One ring per iteration, coupled (processor, cyclic) patches synced
    // every round.
    const label nTotal = nLayers_ + nRampLayers_;
    if (nTotal > 0)
    {
        const labelUList& own = mesh.faceOwner();
        const labelUList& nei = mesh.faceNeighbour();

        DynamicList<label> seedFaces(256);

        for (label facei = 0; facei < mesh.nInternalFaces(); ++facei)
        {
            if (isMixed.test(own[facei]) || isMixed.test(nei[facei]))
            {
                seedFaces.append(facei);
            }
        }
        for
        (
            label facei = mesh.nInternalFaces();
            facei < mesh.nFaces();
            ++facei
        )
        {
            if (isMixed.test(own[facei]))
            {
                seedFaces.append(facei);
            }
        }

        const List<topoDistanceData<label>> seedInfo
        (
            seedFaces.size(),
            topoDistanceData<label>(1, 0)
        );

        List<topoDistanceData<label>> faceDist(mesh.nFaces());

        FaceCellWave<topoDistanceData<label>> wave(mesh, faceDist, cellDist);
        wave.setFaceInfo(seedFaces, seedInfo);
        wave.iterate(nTotal);
    }

    // Hard core: band = 1. Ramp rings: cosine decay towards the bulk.
    // Ramp cells are pure-phase, so fractional weights there stay
    // mass-momentum consistent.
    const scalar pi = constant::mathematical::pi;

    forAll(band_, celli)
    {
        const label d = cellDist[celli].distance();

        if (isMixed.test(celli) || (d >= 1 && d <= nLayers_))
        {
            band_[celli] = 1.0;
        }
        else if (d > nLayers_ && d <= nTotal)
        {
            band_[celli] =
                0.5*(1.0 + Foam::cos(pi*(d - nLayers_)/(nRampLayers_ + 1)));
        }
        else
        {
            band_[celli] = 0.0;
        }
    }
    band_.correctBoundaryConditions();

    // Face blend factor: the larger of the two adjacent cell values, so
    // every face of a core cell is fully upwind.
    {
        const labelUList& own = mesh.faceOwner();
        const labelUList& nei = mesh.faceNeighbour();

        for (label facei = 0; facei < mesh.nInternalFaces(); ++facei)
        {
            faceMask_[facei] = max(band_[own[facei]], band_[nei[facei]]);
        }

        surfaceScalarField::Boundary& bMask = faceMask_.boundaryFieldRef();
        const volScalarField::Boundary& bBand = band_.boundaryField();

        forAll(bMask, patchi)
        {
            scalarField& pMask = bMask[patchi];
            const scalarField ownVals(bBand[patchi].patchInternalField());

            if (mesh.boundary()[patchi].coupled())
            {
                const scalarField neiVals(bBand[patchi].patchNeighbourField());
                forAll(pMask, facei)
                {
                    pMask[facei] = max(ownVals[facei], neiVals[facei]);
                }
            }
            else
            {
                pMask = ownVals;
            }
        }
    }
}


void Foam::interfaceBand::updateIfNeeded()
{
    const label curTimeIndex = band_.mesh().time().timeIndex();

    if (timeIndex_ != curTimeIndex)
    {
        rebuild();
        timeIndex_ = curTimeIndex;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::interfaceBand::interfaceBand(const fvMesh& mesh)
:
    MeshObject<fvMesh, Foam::TopologicalMeshObject, interfaceBand>(mesh),
    timeIndex_(-1),
    nLayers_(-1),
    nRampLayers_(-1),
    threshold_(-1),
    alphaName_(),
    band_
    (
        IOobject
        (
            "interfaceBandCells",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dimless, Zero)
    ),
    faceMask_
    (
        IOobject
        (
            "interfaceBandFaceMask",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        mesh,
        dimensionedScalar(dimless, Zero)
    )
{
    readControls();
}


// * * * * * * * * * * * * * * * * Selector  * * * * * * * * * * * * * * * * //

Foam::interfaceBand& Foam::interfaceBand::New(const fvMesh& mesh)
{
    auto* ptr = mesh.thisDb().getObjectPtr<interfaceBand>("interfaceBand");

    if (!ptr)
    {
        ptr = new interfaceBand(mesh);
        regIOobject::store(ptr);
    }

    return *ptr;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::volScalarField& Foam::interfaceBand::band()
{
    updateIfNeeded();
    return band_;
}


const Foam::surfaceScalarField& Foam::interfaceBand::faceMask()
{
    updateIfNeeded();
    return faceMask_;
}


// ************************************************************************* //
