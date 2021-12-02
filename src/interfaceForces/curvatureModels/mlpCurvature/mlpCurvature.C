/*---------------------------------------------------------------------------*\
            Copyright (c) 2021, Asim Onder
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

#include "mlpCurvature.H"
#include "addToRunTimeSelectionTable.H"

//#include "alphaContactAngleFvPatchScalarField.H"
#include "alphaContactAngleTwoPhaseFvPatchScalarField.H"
#include "mathematicalConstants.H"
#include "surfaceInterpolate.H"
#include "fvcDiv.H"
#include "fvcGrad.H"
#include "reconstructionSchemes.H"
#include "OFstream.H"
#include "IFstream.H"
#include "unitConversion.H"
#include "cellSet.H"
#include <math.h>
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(mlpCurvature, 0);
    addToRunTimeSelectionTable(curvatureModel,mlpCurvature, components);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mlpCurvature::mlpCurvature
(
    const dictionary& dict,
    const volScalarField& alpha1,
    const volVectorField& U
)
:
    curvatureModel
    (
        typeName,
	dict,
        alpha1,
        U
    ),
    deltaN_
    (
        "deltaN",
        1e-8/pow(average(alpha1.mesh().V()), 1.0/3.0)
     ),
    mesh_(alpha1.mesh()),
    RDF_(alpha1.mesh()),
    //ijkMesh_(alpha1.mesh()),
    ijkMesh_(mesh_),
    globalNumbering_(ijkMesh_.globalNumbering()),
    boundaryCells_(mesh_.nCells(),false),
    //bias_(dict.lookupOrDefault<bool>("use_bias",false)),
    zonalModel_(dict.lookupOrDefault<bool>("zonalModel",true)),
    iMax_(dict.lookupOrDefault<label>("iMax",1)),
    jMax_(dict.lookupOrDefault<label>("jMax",1)),
    kMax_(dict.lookupOrDefault<label>("kMax",2)),
    //stencil_(mesh_,ijkMesh_,2,2,1)
    stencil_(mesh_,ijkMesh_,2)
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

  //is2D_=false;
  //if(ijkMesh_.isEmpty().x() or ijkMesh_.isEmpty().y() or ijkMesh_.isEmpty().z())
  Info<<"ijkMesh_.isEmpty()="<<ijkMesh_.isEmpty()<<endl;
  is2D_=true;

  if(!ijkMesh_.isEmpty().x())
    dl_=dx;
  else
    {
      if(!ijkMesh_.isEmpty().y())
	dl_=dy;
      else
	dl_=dz;
    }

  Info<<"mlpCurvature: grid spacing, dl_="<<dl_<<endl;

  nMax_=max(iMax_,kMax_);
  if (!is2D_)
    nMax_=max(nMax_,jMax_);
  
  ijkMesh_.markBoundaryCells(boundaryCells_,nMax_);
  Info<<"mlpCurvature: done setting up ijkMesh..."<<endl;

  std::string fName=mesh_.time().path()/"machineLearningModels/curv";
  if (Pstream::parRun())
    fName=mesh_.time().path()/"../machineLearningModels/curv";

  mlp_=Foam::multilayerPerceptron::multilayerPerceptron(fName);
  NInput_=mlp_.stencilSize();
  Info<<"mlpCurvature:stencil size="<<NInput_<<endl;
  if (is2D_)
    {
      if (NInput_!=(2*iMax_+1)*(2*kMax_+1))
	{
	  Info<<"iMax="<<iMax_<<"; kMax_="<<kMax_<<"; input layer="<<NInput_<<endl;
	  FatalErrorInFunction
	    << "Stencil size is inconsistent with the mlp input layer dimension!"
	    << abort(FatalError);
	}
    }
  else
    {
      if (NInput_!=(2*iMax_+1)*(2*jMax_+1)*(2*kMax_+1))
	{
	  Info<<"iMax="<<iMax_<<"; jMax_="<<jMax_<<"; kMax_="<<kMax_<<"; input layer="<<NInput_<<endl;
	  FatalErrorInFunction
	    << "Stencil size is inconsistent with mlp input layer dimension!"
	    << abort(FatalError);
	}
    }

    

  //  std::vector<int>
  indices_.resize(NInput_);
  int iA=0;
  if (is2D_)
    {      
      //nMax=max(jMax,kMax);
      //hard coded to z as 2nd direction
      Info<<"mlpCurvature currently supports only the x-z coordinate system in 2D."<<endl;
      for (int i=-iMax_;i<iMax_+1;i++)
	{
	  for (int k=-kMax_;k<kMax_+1;k++)
	    {
	      indices_[iA]=stencil_.a2(i,k);
	      iA+=1;
	    }
	}
    }
  else
    {
      for (int i=-iMax_;i<iMax_+1;i++)
	{
	  for (int j=-jMax_;j<jMax_+1;j++)
	    {
	      for (int k=-kMax_;k<kMax_+1;k++)
		{
		  indices_[iA]=stencil_.a3(i,j,k);
		  iA+=1;
		}
	    }
	}
    }
}

// * * * * * * * * * * * * * * Public Access Member Functions  * * * * * * * * * * * * * * //
void Foam::mlpCurvature::correctContactAngle
(
    surfaceVectorField::Boundary& nHatb,
    const surfaceVectorField::Boundary& gradAlphaf
) const
{
    const fvMesh& mesh = alpha1_.mesh();
    const volScalarField::Boundary& abf = alpha1_.boundaryField();

    const fvBoundaryMesh& boundary = mesh.boundary();

    forAll(boundary, patchi)
    {
        if (isA<alphaContactAngleTwoPhaseFvPatchScalarField>(abf[patchi]))
        {
            alphaContactAngleTwoPhaseFvPatchScalarField& acap =
                const_cast<alphaContactAngleTwoPhaseFvPatchScalarField&>
                (
                    refCast<const alphaContactAngleTwoPhaseFvPatchScalarField>
                    (
                        abf[patchi]
                    )
                );

            fvsPatchVectorField& nHatp = nHatb[patchi];
            const scalarField theta
            (
                degToRad() * acap.theta(U_.boundaryField()[patchi], nHatp)
            );

            const vectorField nf
            (
                boundary[patchi].nf()
            );

            // Reset nHatp to correspond to the contact angle

            const scalarField a12(nHatp & nf);
            const scalarField b1(cos(theta));

            scalarField b2(nHatp.size());
            forAll(b2, facei)
            {
                b2[facei] = cos(acos(a12[facei]) - theta[facei]);
            }

            const scalarField det(1.0 - a12*a12);

            scalarField a((b1 - a12*b2)/det);
            scalarField b((b2 - a12*b1)/det);

            nHatp = a*nf + b*nHatp;
            nHatp /= (mag(nHatp) + deltaN_.value());

            acap.gradient() = (nf & nHatp)*mag(gradAlphaf[patchi]);
            acap.evaluate();
        }
    }
}

Foam::scalar Foam::mlpCurvature::estimateSignK(const List<scalar>& A, const vector& n)
{
  bool posK=false;
  if (is2D_)
    {      
      List<scalar> H;
      H.setSize(3,0.0);
      H=0;
      label dir=2;
      if (Foam::mag(n.x())>Foam::mag(n.z()))
	dir=0;
      
      if (dir==2)
	{
	  for (int i=-1;i<2;i++)
	    {
	      for (int k=-kMax_;k<kMax_+1;k++)
		{
		  H[i+1]+=A[stencil_.a2(i,k)];
		}
	    }
	}
      else if(dir==0)
	{
	  for (int k=-1;k<2;k++)
	    {
	      for (int i=-iMax_;i<iMax_+1;i++)
		{
		  H[k+1]+=A[stencil_.a2(i,k)];
		}
	    }
	}
      
      scalar Hxx=(H[0]-2.*H[1]+H[2]);
      posK=!std::signbit(-Hxx);
    }
  else
    NotImplemented;

  return (posK?1.0:-1.0);

}

void Foam::mlpCurvature::calculateK()
{
  Info<<"mlpCurvature: calculating the curvature..."<<endl;
  const fvMesh& mesh = alpha1_.mesh();
  const volVectorField& C = mesh.C();
  //const surfaceVectorField& Sf = mesh.Sf();
  point Pmin=ijkMesh_.Pmin();
  point Pmax=ijkMesh_.Pmax();
  const scalar dx=ijkMesh_.dx();
  const scalar dy=ijkMesh_.dy();
  const scalar dz=ijkMesh_.dz();
  const label Nx=ijkMesh_.Nx();
  const label Ny=ijkMesh_.Ny();
  const label Nz=ijkMesh_.Nz();
  const labelList& globalIds=ijkMesh_.globalIds(); 
  //modify to dynamic version for speed up
  Map<scalar> alphaIJK;
  Map<vector> faceCentreIJK;
  Map<scalar> curvIJK;

  // where to put this one?
  const volVectorField gradAlpha(fvc::grad(alpha1_, "nHat"));
  surfaceVectorField gradAlphaf(fvc::interpolate(gradAlpha));
  surfaceVectorField nHatfv(gradAlphaf/(mag(gradAlphaf) + deltaN_));
  correctContactAngle(nHatfv.boundaryFieldRef(), gradAlphaf.boundaryFieldRef());

  forAll(K_,celli)
    K_[celli]=0;

  reconstructionSchemes& surf =
    mesh.lookupObjectRef<reconstructionSchemes>("reconstructionScheme");

  //surf.reconstruct(false); <==== this is causing additional surface reconstruction
  //surf.reconstruct(true);

  const boolList& interfaceCells = surf.interfaceCell();
    
  const volVectorField& faceCentre = surf.centre();
  const volVectorField& faceNormal = surf.normal();

  Vector<label> stencilSize(nMax_,nMax_,nMax_);
  ijkMesh_.getZoneField(interfaceCells,alphaIJK,alpha1_,stencilSize);
  label nBoundaryCells=0;
  label nInterfaceCells=0;
  forAll(interfaceCells,celli)
    {    
      vector n = faceNormal[celli];
      if(interfaceCells[celli] && mag(n)>0)
	{
	  nInterfaceCells+=1;
	  if (boundaryCells_[celli])
	    {
	      //printf("Cell at the boundary: proc=%d, x=%f, y=%f, z=%f\n",Pstream::myProcNo(),C[celli].x(),C[celli].y(),C[celli].z());
	      nBoundaryCells+=1;
	      continue;
	    }

	  std::vector<double> inputData(NInput_,0.0);
	  stencil_.setStencil(alphaIJK,ijkMesh_.ijk3(celli));
	  List<scalar> A=stencil_.getStencil();
	  int iA=0;
	  scalar sgnK=1.0;
	  if (is2D_)
	    {      
	      //hard coded using z as 2nd direction
	      if (zonalModel_)
		{
		  //sgnK=estimateSignK(A,n);
		  sgnK=stencil_.estimateSignK();
		  int iF=1;int kF=1;
		  if (-n.x()<0)
		    iF=-1;
		  if (-n.z()<0)
		    kF=-1;
		  if (sgnK<0)
		    {
		      A=1.-A;
		      iF*=-1;kF*=-1;
		    }

		  if (Foam::mag(n.x())<=Foam::mag(n.z()))
		    {
		      for (int i=-iMax_;i<iMax_+1;i++)
			{
			  for (int k=-kMax_;k<kMax_+1;k++)
			    {
			      indices_[iA]=stencil_.a2(iF*i,kF*k);
			      iA+=1;
			    }
			}
		    }
		  else 
		    {
		      for (int i=-iMax_;i<iMax_+1;i++)
			{
			  for (int k=-kMax_;k<kMax_+1;k++)
			    {
			      indices_[iA]=stencil_.a2(iF*k,kF*i);
			      iA+=1;
			    }
			}
		    }
		}
	      else
		{
		  for (int i=-iMax_;i<iMax_+1;i++)
		    {
		      for (int k=-kMax_;k<kMax_+1;k++)
			{
			  indices_[iA]=stencil_.a2(i,k);
			  iA+=1;
			}
		    }
		}
	      
	      iA=0;
	      for(int i: indices_)
		{
		  inputData[iA]=A[i];
		  iA+=1;
		}
	      K_[celli]=sgnK*mlp_.predict(inputData)[0]/double(dl_);
	    }
	  else
	    {
	      NotImplemented;
	    }
	}      
    }
  if(nBoundaryCells>0)
    printf("%d of %d interface cells were at the boundary in proc=%d\n",nBoundaryCells,nInterfaceCells,Pstream::myProcNo());
  K_.correctBoundaryConditions();
  label neiMax=1;
  //this one can cause problems at cyclic BC!
  RDF_.correctBoundaryConditions();
  RDF_.markCellsNearSurf(interfaceCells,neiMax);
  const boolList& nextToInterface =RDF_.nextToInterface(); 
  
  boolList neiInterface(mesh.nCells(),false);

  forAll(neiInterface,celli)
    {
      if (nextToInterface[celli] and !interfaceCells[celli] and !boundaryCells_[celli])
	  neiInterface[celli]=true;
    }

  neiMax+=1;
  stencilSize=Vector<label>(neiMax,neiMax,neiMax);
  ijkMesh_.getZoneField(neiInterface,faceCentreIJK,faceCentre,stencilSize);
  ijkMesh_.getZoneField(neiInterface,curvIJK,K_,stencilSize);


  label iMax,jMax,kMax;

  if (ijkMesh_.isEmpty().x())
    {
      iMax=0;jMax=neiMax;kMax=neiMax;
    }
  else if (ijkMesh_.isEmpty().y())
    {
      iMax=neiMax;jMax=0;kMax=neiMax;
    }
  else if (ijkMesh_.isEmpty().z())
    {
      iMax=neiMax;jMax=neiMax;kMax=0;
    }
  else
    {
      iMax=neiMax;jMax=neiMax;kMax=neiMax;
    }
    
  //label nCount=0;  
  forAll(neiInterface,celli)
    {
      if (neiInterface[celli])
	{
	  const point cc = C[celli];

	  label i=round((cc.x()-Pmin.x())/dl_);
	  label j=round((cc.y()-Pmin.y())/dl_);
	  label k=round((cc.z()-Pmin.z())/dl_);

	  label gblIdMin=-1;
	  scalar distMin=GREAT;
	  for (int ii=-iMax;ii<iMax+1;ii++)
	    {
	      for (int jj=-jMax;jj<jMax+1;jj++)
		{
		  for (int kk=-kMax;kk<kMax+1;kk++)
		    {
		      label il=(i+ii+Nx)%Nx;
		      label jl=(j+jj+Ny)%Ny;
		      label kl=(k+kk+Nz)%Nz; 
		      label ijk=il+Nx*jl+Nx*Ny*kl;
		      label gblId=globalIds[ijk];
		      vector faceCentre=faceCentreIJK[gblId];
		      scalar dist = mag(cc-faceCentre);
		      if (dist < distMin && mag(faceCentre)>0)
			{
			  distMin = dist;
			  gblIdMin = gblId;
			}
		      
		    }
		}
	    }
	  
	  K_[celli]=curvIJK[gblIdMin];
	}
    }
  K_.correctBoundaryConditions();
}



// ************************************************************************* //
