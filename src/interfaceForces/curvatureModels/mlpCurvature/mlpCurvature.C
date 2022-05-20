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
    ijkMesh_(mesh_),
    globalNumbering_(ijkMesh_.globalNumbering()),
    boundaryCells_(mesh_.nCells(),false),
    zonalModel_(dict.lookupOrDefault<bool>("zonalModel",true)),
    useScaling_(dict.lookupOrDefault<bool>("useScaling",true)),
    fillNeighbours_(dict.lookupOrDefault<label>("fillNeighbours",-1)),
    averageZones_(dict.lookupOrDefault<bool>("averageZones",true)),
    //extendInterface_(dict.lookupOrDefault<bool>("extendInterface",false)),
    mlpModel_(dict.lookupOrDefault<word>("mlpModel","machineLearningModels")),
    interfaceTol_(dict.lookupOrDefault<scalar>("interfaceTol",1e-6)),
    xoffsetInput_(dict.lookupOrDefault<scalar>("xoffsetInput",0)),
    yminInput_(dict.lookupOrDefault<scalar>("yminInput",0)),
    gainInput_(dict.lookupOrDefault<scalar>("gainInput",1.)),
    xoffsetOutput_(dict.lookupOrDefault<scalar>("xoffsetOutput",0)),
    yminOutput_(dict.lookupOrDefault<scalar>("yminOutput",0)),
    gainOutput_(dict.lookupOrDefault<scalar>("gainOutput",1.)),
    iMax_(dict.lookupOrDefault<label>("iMax",1)),
    jMax_(dict.lookupOrDefault<label>("jMax",1)),
    kMax_(dict.lookupOrDefault<label>("kMax",1)),
    nMax_(dict.lookupOrDefault<label>("nMax",2)),
    //stencil_(mesh_,ijkMesh_,2,2,1)
    stencil_(mesh_,ijkMesh_,nMax_) //Fix this later! Cannot instantiate the object!!!!
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

  nMax_=max(iMax_,jMax_);
  nMax_=max(nMax_,kMax_);

  //Foam::uniformStencil stencil_(mesh_,ijkMesh_,nMax_);
  //stencil_=Foam::uniformStencil::uniformStencil(mesh_,ijkMesh_,nMax_);
  
  ijkMesh_.markBoundaryCells(boundaryCells_,nMax_);
  Info<<"mlpCurvature: done setting up ijkMesh..."<<endl;

  std::string fName=mesh_.time().path()/mlpModel_+"/curv";
  if (Pstream::parRun())
    fName=mesh_.time().path()+"/../"+mlpModel_+"/curv";

  mlp_=multilayerPerceptron(fName);
  NInput_=mlp_.stencilSize();
  Info<<"mlpCurvature:stencil size="<<NInput_<<endl;
  if (is2D_)
    {
      if (NInput_!=(2*iMax_+1)*(2*jMax_+1))
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
  //indices_.resize(NInput_);

  /*int iA=0;
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
	}*/
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
  const Vector<bool>& isEmpty=ijkMesh_.isEmpty();
  const labelList& globalIds=ijkMesh_.globalIds();
  const globalIndex& globalNumbering=ijkMesh_.globalNumbering();
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
      scalar K1=0.;
      scalar K2=0.;
      scalar K3=0.;
      scalar K4=0.;
      if(interfaceCells[celli]&& mag(n)>0)// and alpha1_[celli]>interfaceTol_ and alpha1_[celli]<1-interfaceTol_) // && mag(n)>0)
	{
	  nInterfaceCells+=1;
	  if (boundaryCells_[celli])
	    {
	      //printf("Cell at the boundary: proc=%d, x=%f, y=%f, z=%f\n",Pstream::myProcNo(),C[celli].x(),C[celli].y(),C[celli].z());
	      nBoundaryCells+=1;
	      continue;
	    }
	  std::vector<int> indices;
	  std::vector<int> indices2;
	  std::vector<int> indices3;
	  std::vector<int> indices4;
	  std::vector<double> inputData(NInput_,0.0);
	  stencil_.setStencil(alphaIJK,ijkMesh_.ijk3(celli));
	  List<scalar> A=stencil_.getStencil();
	  int iA=0;
	  if (is2D_)
	    {
	      if (zonalModel_)	     
		{
		  A=gainInput_*(A-xoffsetInput_)+yminInput_;
		  //vector n = stencil_.calcYoungNormal();
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
		  if (Foam::mag(nT)<=Foam::mag(nN))
		  //if (Foam::mag(nT)<Foam::mag(nN))
		    {
		      for (int t=-tMax;t<tMax+1;t++)
			{
			  for (int n=-nMax;n<nMax+1;n++)
			    {
			      indices.push_back(stencil_.a2(tF*t,nF*n));
			      if (averageZones_)
				{
				  indices2.push_back(stencil_.a2(-tF*t,nF*n));
				  indices3.push_back(stencil_.a2(tF*t,-nF*n));
				  indices4.push_back(stencil_.a2(-tF*t,-nF*n));
				}
			    }
			}
		    }
		  else 
		    {
		      for (int t=-tMax;t<tMax+1;t++)
			{
			  for (int n=-nMax;n<nMax+1;n++)
			    {
			      indices.push_back(stencil_.a2(tF*n,nF*t));
			      if (averageZones_)
				{
				  indices2.push_back(stencil_.a2(-tF*n,nF*t));
				  indices3.push_back(stencil_.a2(tF*n,-nF*t));
				  indices4.push_back(stencil_.a2(-tF*n,-nF*t));
				}
			    }
			}
		    }

		  for(int i=0;i<NInput_;i++)
		    {
		      inputData[i]=A[indices[i]];
		    }
		  if (averageZones_)
		    {
		      K1=mlp_.predict(inputData)[0];
		      K_[celli]=0.25*K1; //mlp_.predict(inputData)[0];///double(dl_);
		      for(int i=0;i<NInput_;i++)
			{
			  inputData[i]=A[indices2[i]];
			}
		      K2=mlp_.predict(inputData)[0];
		      K_[celli]+=0.25*K2;
		      for(int i=0;i<NInput_;i++)
			{
			  inputData[i]=A[indices3[i]];
			}
		      K3=mlp_.predict(inputData)[0];
		      K_[celli]+=0.25*K3;
		      for(int i=0;i<NInput_;i++)
			{
			  inputData[i]=A[indices4[i]];
			}
		      K4=mlp_.predict(inputData)[0];
		      K_[celli]+=0.25*K4;
		    }
		  else
		    {
		      K_[celli]=mlp_.predict(inputData)[0];
		    }

		  K_[celli]=(K_[celli]-yminOutput_)/gainOutput_+xoffsetOutput_;
		  K_[celli]/=double(dl_);
		  indices.clear();
		  if (averageZones_)
		    {
		      indices2.clear();
		      indices3.clear();
		      indices4.clear();
		    }
		}
	      else
		{
		  if (useScaling_)
		    {
		      A=gainInput_*(A-xoffsetInput_)+yminInput_;
		    }
		  for (int i=-iMax_;i<iMax_+1;i++)
		    {
		      for (int j=-jMax_;j<jMax_+1;j++)
			{
			  indices.push_back(stencil_.a2(i,j));
			}
		    }
		  for(int i=0;i<NInput_;i++)
			{
			  inputData[i]=A[indices[i]];
			}
		  K_[celli]=mlp_.predict(inputData)[0];
		  if(useScaling_)
		    {
		      K_[celli]=(K_[celli]-yminOutput_)/gainOutput_+xoffsetOutput_;
		    }
		  K_[celli]/=double(dl_);

		}
	    }
	  else
	    {
	      NotImplemented;
	    }
	}      
    }
  if(nBoundaryCells>0)
    printf("%d of %d interface cells were at the boundary in proc=%d\n",nBoundaryCells,nInterfaceCells,Pstream::myProcNo());


  #include "../interpolateCutCellsToFaces.H"
  
      
}



// ************************************************************************* //
