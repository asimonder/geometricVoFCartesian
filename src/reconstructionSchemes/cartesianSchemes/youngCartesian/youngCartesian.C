/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2021 Asim Onder
-------------------------------------------------------------------------------
License
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

#include "youngCartesian.H"
#include "fvc.H"
#include "addToRunTimeSelectionTable.H"
#include "globalIndex.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace reconstruction
{
    defineTypeNameAndDebug(youngCartesian, 0);
    addToRunTimeSelectionTable(reconstructionSchemes, youngCartesian, components);
}
}



void Foam::reconstruction::youngCartesian::gradSurf(const volScalarField& phi)
{

  label iMax=1;
  Map<scalar> phiIJK;
  Vector<label> stencilSize(iMax,iMax,iMax);
  ijkMesh_.getZoneField(interfaceCell_,phiIJK,phi,stencilSize);

  forAll(interfaceLabels_, i)
    {
      label celli=interfaceLabels_[i];
      stencil_.setStencil(phiIJK,ijkMesh_.ijk3(celli));
      interfaceNormal_[i] = stencil_.calcYoungNormal();      
    }
  
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::reconstruction::youngCartesian::youngCartesian
(
    volScalarField& alpha1,
    const surfaceScalarField& phi,
    const volVectorField& U,
    const dictionary& dict
)
:
    reconstructionSchemes
    (
        typeName,
        alpha1,
        phi,
        U,
        dict
    ),
    mesh_(alpha1.mesh()),
    interfaceNormal_(fvc::grad(alpha1)),
    isoFaceTol_(modelDict().lookupOrDefault<scalar>("isoFaceTol", 1e-8)),
    surfCellTol_(modelDict().lookupOrDefault<scalar>("surfCellTol", 1e-8)),
    sIterPLIC_(mesh_,surfCellTol_),
    ijkMesh_(mesh_),
    stencil_(mesh_,ijkMesh_,1),
    boundaryCells_(mesh_.nCells(),false)
{
  ijkMesh_.markBoundaryCells(boundaryCells_,1);
  reconstruct();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::reconstruction::youngCartesian::reconstruct(bool forceUpdate)
{
    const bool uptodate = alreadyReconstructed(forceUpdate);

    if (uptodate && !forceUpdate)
    {
        return;
    }

    if (mesh_.topoChanging())
    {
        // Introduced resizing to cope with changing meshes
        //if(interfaceCell_.size() != mesh_.nCells())
        //{
        //    interfaceCell_.resize(mesh_.nCells());
        //}
      Info<<"Dynamic mesh is not supported!"<<endl;
      Foam::FatalError();
    }
    
    interfaceCell_ = false;

    interfaceLabels_.clear();

    forAll(alpha1_,celli)
    {
        if(sIterPLIC_.isASurfaceCell(alpha1_[celli]))
        {
            interfaceCell_[celli] = true; // is set to false earlier
            interfaceLabels_.append(celli);
        }
    }
    interfaceNormal_.setSize(interfaceLabels_.size());
    centre_ = dimensionedVector("centre", dimLength, vector::zero);
    normal_ = dimensionedVector("normal", dimArea, vector::zero);

    gradSurf(alpha1_);

    forAll(interfaceLabels_, i)
    {
        const label celli = interfaceLabels_[i];
        if (mag(interfaceNormal_[i]) == 0)
        {
            continue;
        }

        sIterPLIC_.vofCutCell
        (
            celli,
            alpha1_[celli],
            isoFaceTol_,
            100,
            interfaceNormal_[i]
        );

        if (sIterPLIC_.cellStatus() == 0)
        {
            normal_[celli] = sIterPLIC_.surfaceArea();
            centre_[celli] = sIterPLIC_.surfaceCentre();
            if (mag(normal_[celli]) == 0)
            {
                normal_[celli] = vector::zero;
                centre_[celli] = vector::zero;
            }

        }
        else
        {
            normal_[celli] = vector::zero;
            centre_[celli] = vector::zero;
        }
    }
}

void Foam::reconstruction::youngCartesian::mapAlphaField() const
{
    // without it, we seem to get a race condition
    mesh_.C();

    cutCellPLIC cutCell(mesh_);

    forAll(normal_, celli)
    {
        if (mag(normal_[celli]) != 0)
        {
            vector n = normal_[celli]/mag(normal_[celli]);
            scalar cutValue = (centre_[celli] - mesh_.C()[celli]) & (n);
            cutCell.calcSubCell
            (
                celli,
                cutValue,
                n
            );
            alpha1_[celli] = cutCell.VolumeOfFluid();

        }
    }
    alpha1_.correctBoundaryConditions();
    alpha1_.oldTime () = alpha1_;
    alpha1_.oldTime().correctBoundaryConditions();

}


// ************************************************************************* //
