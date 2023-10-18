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

#include "mlpNormal.H"
#include "fvc.H"
#include "addToRunTimeSelectionTable.H"
#include "globalIndex.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace reconstruction
{
    defineTypeNameAndDebug(mlpNormal, 0);
    addToRunTimeSelectionTable(reconstructionSchemes, mlpNormal, components);
}
}


void Foam::reconstruction::mlpNormal::gradSurf(const volScalarField& phi)
{
  Map<scalar> phiIJK;
  Vector<label> stencilSize(iMax_,jMax_,kMax_);
  ijkMesh_.getZoneField(interfaceCell_,phiIJK,phi,stencilSize);
  std::vector<double> inputData(NInput_,0.0);
  const Vector<bool>& isEmpty=ijkMesh_.isEmpty();
  
  //Info<<"iMax_="<<iMax_<<", kMax_="<<kMax_<<", nMax_="<<kMax_<<endl;
  
  /*Info<<"iMax_="<<iMax_<<", kMax_="<<kMax_<<", nMax_="<<kMax_<<endl;
  for (int i=-1;i<1+1;i++)
    {
      for (int k=-1;k<1+1;k++)
	{
	  Info<<stencil_.a2(i,k)<<endl;
	}
	}*/
  forAll(interfaceLabels_, iCell)
    {
      //Info<<"iCell="<<iCell<<endl;
      vector m(0.0,0.0,0.0);
      label celli=interfaceLabels_[iCell];
      stencil_.setStencil(phiIJK,ijkMesh_.ijk3(celli));
      List<scalar> A=stencil_.getStencil();
      int iA=0;	      
      scalar sgnK=1.0;
      if (is2D_)
	{      
	  if (zonalModel_)
	    {
	      vector n = stencil_.calcYoungNormal();
	      sgnK=0; //stencil_.estimateSignK(n);
	      scalar nN=0;
	      scalar nT=0;
	      label tMax=0;
	      label nMax=0;
	      if (isEmpty.x())
		{
		  nT=n.y();
		  nN=n.z();
		}
	      else if (isEmpty.y())
		{
		  nT=n.x();
		  nN=n.z();
		}
	      else 
		{
		  nT=n.x();
		  nN=n.y();
		}
	      tMax=iMax_;
	      nMax=jMax_;
	      
	      int tF=1;int nF=1;
	      if (-nT<0)
		tF=-1;
	      if (-nN<0)
		nF=-1;
	      if (sgnK<0)
		{
		  A=1.-A;
		  tF*=-1;nF*=-1;
		}
	      
 
	      if (Foam::mag(nT)<=Foam::mag(nN))
		{
		  for (int t=-tMax;t<tMax+1;t++)
		    {
		      for (int n=-nMax;n<nMax+1;n++)
			{
			  indices_[iA]=stencil_.a2(tF*t,nF*n);
			  iA+=1;
			}
		    }
		}
	      else 
		{
		  for (int t=-tMax;t<tMax+1;t++)
		    {
		      for (int n=-nMax;n<nMax+1;n++)
			{
			  indices_[iA]=stencil_.a2(tF*n,nF*t);
			  iA+=1;
			}
		    }
		}

	      
	      iA=0;
	      for(int ii: indices_)
		{
		  inputData[iA]=A[ii];
		  iA+=1;
		}

	      scalar nx=mlp_.predict(inputData)[0];
	      scalar nMag=Foam::sqrt(Foam::sqr(nx)+1);
	      m.x()=nx/nMag;
	      m.z()=1./nMag;
	      if (Foam::mag(n.x())<=Foam::mag(n.z()))
		{		  
		  m.x()*=std::copysign(1.0,n.x()); 
		  m.z()*=std::copysign(1.0,n.z()); 
		}
	      else
		{
		  vector mR=m;
		  m.x()=std::copysign(1.0,n.x())*mR.z(); 
		  m.z()=std::copysign(1.0,n.z())*mR.x(); 
		}		  
	    }
	  else
	    {
	      iA=0;
	      for (int i=-iMax_;i<iMax_+1;i++)
		{
		  for (int k=-kMax_;k<kMax_+1;k++)
		    {
		      indices_[iA]=stencil_.a2(i,k);
		      iA+=1;
		    }
		}
	      iA=0;
	      for(int ii: indices_)
		{
		  inputData[iA]=A[ii];
		  iA+=1;
		}
	      scalar nx=mlp_.predict(inputData)[0];
	      scalar nMag=Foam::sqrt(Foam::sqr(nx)+1);
	      m.x()=nx/nMag;
	      m.z()=1./nMag;	      
	    }
	  interfaceNormal_[iCell] = m;       	      
	}
      else
	NotImplemented;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
Foam::reconstruction::mlpNormal::mlpNormal
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
    boundaryCells_(mesh_.nCells(),false),
    //bias_(modelDict().lookupOrDefault<bool>("use_bias",false)),
    zonalModel_(dict.lookupOrDefault<bool>("zonalModel",true)),
    iMax_(dict.lookupOrDefault<label>("iMax",1)),
    jMax_(dict.lookupOrDefault<label>("jMax",1)),
    kMax_(dict.lookupOrDefault<label>("kMax",1)),
    nMax_(1),
    stencil_(mesh_,ijkMesh_,1)
{
  const scalar dx=ijkMesh_.dx();
  const scalar dy=ijkMesh_.dy();
  const scalar dz=ijkMesh_.dz();

  bool isUniform=true;
  if ( !ijkMesh_.isEmpty().x() and !ijkMesh_.isEmpty().y() and fabs(1-dx/dy)/dx>1e-6)
    isUniform=false;
  if ( !ijkMesh_.isEmpty().y() and !ijkMesh_.isEmpty().z() and fabs(1-dy/dz)/dy>1e-6)
    isUniform=false;
  if ( !ijkMesh_.isEmpty().x() and !ijkMesh_.isEmpty().z() and fabs(1-dx/dz)/dx>1e-6)
    isUniform=false;

  if (!isUniform)
    FatalErrorInFunction
      << "mlpCurvature requires uniform grids!"
      << abort(FatalError);

  is2D_=false;
  if(ijkMesh_.isEmpty().x() or ijkMesh_.isEmpty().y() or ijkMesh_.isEmpty().z())
    is2D_=true;

  //nMax_=max(iMax_,kMax_);
  //if (!is2D_)
  // nMax_=max(nMax_,jMax_);

  //stencil_=Foam::uniformStencil::uniformStencil(mesh_,ijkMesh_,nMax_);

  if (zonalModel_)
    Info<<"zonalModel=true."<<endl;
  ijkMesh_.markBoundaryCells(boundaryCells_,nMax_);

  std::string fName=mesh_.time().path()/"machineLearningModels/der1";
  if (Pstream::parRun())
    fName=mesh_.time().path()/"../machineLearningModels/der1";

  mlp_=multilayerPerceptron(fName);
  NInput_=mlp_.stencilSize();

  is2D_=false;
  if(ijkMesh_.isEmpty().x() or ijkMesh_.isEmpty().y() or ijkMesh_.isEmpty().z())
    is2D_=true;

  Info<<"mlpNormal:stencil size="<<NInput_<<endl;
  if (is2D_)
    {
      if (NInput_!=(2*iMax_+1)*(2*kMax_+1))
	{
	  Info<<"iMax="<<iMax_<<"; kMax_="<<kMax_<<"; input layer="<<NInput_<<endl;
	  FatalErrorInFunction
	    << "Stencil size is inconsistent with mlp input dimension of the mlpNormal model!"
	    << abort(FatalError);
	}
    }
  else
    {
      if (NInput_!=(2*iMax_+1)*(2*jMax_+1)*(2*kMax_+1))
	{
	  FatalErrorInFunction
	    << "Stencil size is inconsistent with mlp input dimension of the mlpNormal!"
	    << abort(FatalError);
	}
    }
  
  indices_.resize(NInput_);
  //int iA=0;

  //if(mlp_.stencilSize()=!(2*iMax_+1)*(2*jMax_+1)*(kMax_+1))
  //  Foam::FatalErrorInFunction("Dimension mismatch between the MLP and specified stencil size.");
    
  
  /*if(mlp_.stencilSize()==27)
    {
      iMax_=1;jMax_=1;kMax_=1;
      stencil_=Foam::uniformStencil(mesh_,ijkMesh_,iMax_);
    }
  else if (mlp_.stencilSize()==75)
    {
      iMax_=2;jMax_=2;kMax_=1;
      stencil_=Foam::uniformStencil(mesh_,ijkMesh_,iMax_,jMax_,kMax_);
    }
  else if (mlp_.stencilSize()==125)
    {
      iMax_=2;jMax_=2;kMax_=2;
      stencil_=Foam::uniformStencil(mesh_,ijkMesh_,iMax_);
    }
  else
  Info<<"stencil size cannot be mapped to 3D!"<<endl;*/

  Info<<"Reconstructing the interface..."<<endl;
      
  reconstruct();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void Foam::reconstruction::mlpNormal::reconstruct(bool forceUpdate)
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

void Foam::reconstruction::mlpNormal::mapAlphaField() const
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
