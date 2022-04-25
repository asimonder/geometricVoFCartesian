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
    This file is part of the unofficial extension library geometricVoFCartesian.

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

#include "mycCartesian.H"
#include "fvc.H"
#include "addToRunTimeSelectionTable.H"
#include "globalIndex.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace reconstruction
{
    defineTypeNameAndDebug(mycCartesian, 0);
    addToRunTimeSelectionTable(reconstructionSchemes, mycCartesian, components);
}
}



void Foam::reconstruction::mycCartesian::gradSurf(const volScalarField& phi)
{
  Map<scalar> phiIJK;
  Vector<label> stencilSize(1,1,1);
  ijkMesh_.getZoneField(interfaceCell_,phiIJK,phi,stencilSize);

  //printf("Proc=%d: Done getZoneField.\n",Pstream::myProcNo());
  label youngCount=0;
  label CCDCount=0;
  //const scalar eps=10e-14;
  forAll(interfaceLabels_, i)
    {      
      int celli=interfaceLabels_[i];
      
      if (boundaryCells_[celli])
	{
	  printf("Proc=%d: Interface is at a non-cyclic cellSet or domain boundary. Non-cyclic boundaries is not supported at the moment. Setting the normal of the interfacial the segment in cell %d to zero!\n",Pstream::myProcNo(),celli);
	  //FatalErrorInFunction
	  //  << "Interface is at a non-cyclic cellSet or domain boundary. Non-cyclic boundaries is not supported at the moment."
	  //<< abort(FatalError);
	  interfaceNormal_[i] =vector::zero;
	  continue;
	  }
      stencil_.setStencil(phiIJK,ijkMesh_.ijk3(celli));
      //interfaceNormal_[i] = stencil_.calcCCDNormal();
      
      vector m1=stencil_.calcYoungNormal();
      vector m2=stencil_.calcCCDNormal();
      vector m2C=m2;

      if(Foam::mag(m1)+Foam::mag(m2)==0)
	{
	  printf("Proc=%d: Young and CDC normals are zero! Skipping the interface cell=%d!\n",Pstream::myProcNo(),celli);
	  //Info<< "Young and CDC normals are zero! Skipping the interface cell...";
	  continue;
	}
      
      if(Foam::mag(m1)>0)
	{
	  m1/=(fabs(m1.x())+fabs(m1.y())+fabs(m1.z()));
	}
      if(Foam::mag(m2)>0)
	{
	  m2/=(fabs(m2.x())+fabs(m2.y())+fabs(m2.z()));
	}

      //printf("Proc=%d: done scaling normals.\n",Pstream::myProcNo());

      if (fabs(m2C.y())==1.0)
	{
	  if (fabs(m1.y())<fabs(m2.y()))
	    {
	      interfaceNormal_[i] = m1;
	      youngCount+=1;
	    }
	  else
	    {
	      interfaceNormal_[i] = m2;
	      CCDCount+=1;
	    }
	}
      else if (fabs(m2C.z())==1.0)
	{
	  if (fabs(m1.z())<fabs(m2.z()))
	    {
	      interfaceNormal_[i] = m1;
	      youngCount+=1;
	    }
	  else
	    {
	      interfaceNormal_[i] = m2;
	      CCDCount+=1;
	    }
	}	  
      else if (fabs(m2C.x())==1.0)
	{
	  if (fabs(m1.x())<fabs(m2.x()))
	    {
	      interfaceNormal_[i] = m1;
	      youngCount+=1;
	    }
	  else
	    {
	      interfaceNormal_[i] = m2;
	      CCDCount+=1;
	    }
	}
      else
	FatalErrorInFunction
	  << "Unity direction cannot be found in CCD scheme!"
	  << abort(FatalError);



    }

  //printf("Proc=%d: Done gradSurf.\n",Pstream::myProcNo());

  
  if (Pstream::parRun())
      {
	Foam::reduce(youngCount, sumOp<label>());
	Foam::reduce(CCDCount, sumOp<label>());    
      }

  Info<<"youngCount="<<youngCount<<endl;
  Info<<"CCDCount="<<CCDCount<<endl;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::reconstruction::mycCartesian::mycCartesian
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
  Info<<"Constructing the interface using mycCartesian method..."<<endl;
  ijkMesh_.markBoundaryCells(boundaryCells_,1);
  reconstruct();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::reconstruction::mycCartesian::reconstruct(bool forceUpdate)
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

void Foam::reconstruction::mycCartesian::mapAlphaField() const
{
    // without this line, we seem to get a race condition
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
