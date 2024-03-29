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

#include "heightFunction.H"
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
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(heightFunction, 0);
    addToRunTimeSelectionTable(curvatureModel,heightFunction, components);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::heightFunction::heightFunction
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
    mesh_(alpha1.mesh()),
    deltaN_
    (
        "deltaN",
        1e-8/pow(average(alpha1.mesh().V()), 1.0/3.0)
     ),
    ijkMesh_(alpha1.mesh()),
    globalNumbering_(ijkMesh_.globalNumbering()),
    nMax_(dict.subDict("heightFunctionParams").lookupOrDefault<label>("nMax",3)),
    //fillNeighbours_(dict.lookupOrDefault<label>("fillNeighbours",-1)),
    interfaceTol_(1e-12),
    //rdfMark_(dict.lookupOrDefault<bool>("rdfMark",false)),
    //extendInterface_(false),
    boundaryCells_(mesh_.nCells(),false)
{
  
  ijkMesh_.markBoundaryCells(boundaryCells_,nMax_);
  Info<<"Done setting up ijkMesh for the height function method..."<<endl;

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
      << "Height function requires uniform grids!"
      << abort(FatalError);

  is2D_=false;
  if(ijkMesh_.isEmpty().x() or ijkMesh_.isEmpty().y() or ijkMesh_.isEmpty().z())
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

  Info<<"grid spacing, dl_="<<dl_<<endl;
}

// * * * * * * * * * * * * * * Public Access Member Functions  * * * * * * * * * * * * * * //
void Foam::heightFunction::correctContactAngle
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


Foam::List<Foam::scalar> Foam::heightFunction::calculateHeights(const Map<scalar>& phiIJK,const Vector<label>& ijk,label dir)
{
  const labelList& globalIds=ijkMesh_.globalIds();
  const Vector<bool>& isEmpty=ijkMesh_.isEmpty();
  const label Nx=ijkMesh_.Nx();
  const label Ny=ijkMesh_.Ny();
  const label Nz=ijkMesh_.Nz();
  
  const label& iP=ijk.x();
  const label& jP=ijk.y();
  const label& kP=ijk.z();

  labelList nMax;
  nMax.setSize(3,1);
  nMax[dir]=nMax_;

  int il,jl,kl=0;

  List<scalar> H;
  if (is2D_)
   {
     H.setSize(3,0.0);
     if (isEmpty.x())
       {
	 il=0;
	 for (int j=-nMax[1];j<nMax[1]+1;j++)
	   {
	     if (jP+j+1>Ny and ijkMesh_.symYOut()) 
	       jl=jP+1-j;
	     else if (jP+j<0 and ijkMesh_.symYIn())
	       jl=jP-1-j;
	     else
	       jl=(jP+j+Ny)%Ny;
	     for (int k=-nMax[2];k<nMax[2]+1;k++)
	       {
		 if (kP+k+1>Nz and ijkMesh_.symZOut()) 
		   kl=kP+1-k;
		 else if (kP-k<0 and ijkMesh_.symZIn())
		   kl=kP-1-k;
		 else
		   kl=(kP+k+Nz)%Nz;
		 label ijkG=ijkMesh_.ijk1(il,jl,kl); 
		 label gblIdx=globalIds[ijkG];
		 if (dir==1)
		   H[k+1]+=phiIJK[gblIdx];
		 else
		   H[j+1]+=phiIJK[gblIdx];
	       }
	   }
       }
     else if (isEmpty.y())
       {
	 jl=0;
	 for (int i=-nMax[0];i<nMax[0]+1;i++)
	   {
	     if (iP+i+1>Nx and ijkMesh_.symXOut())
	       il=iP+1-i;
	     else if (iP+i<0 and ijkMesh_.symXIn())
	       il=iP-1-i;
	     else
	       il=(iP+i+Nx)%Nx;
	     for (int k=-nMax[2];k<nMax[2]+1;k++)
	       {
		 if (kP+k+1>Nz and ijkMesh_.symZOut()) 
		   kl=kP+1-k;
		 else if (kP-k<0 and ijkMesh_.symZIn())
		   kl=kP-1-k;
		 else
		   kl=(kP+k+Nz)%Nz;
		 label ijkG=ijkMesh_.ijk1(il,jl,kl); 
		 label gblIdx=globalIds[ijkG];
		 if (dir==0)
		   H[k+1]+=phiIJK[gblIdx];
		 else
		   H[i+1]+=phiIJK[gblIdx];
	       }
	   }
       }
     else 
       {
	 kl=0;
	 for (int i=-nMax[0];i<nMax[0]+1;i++)
	   {
	     if (iP+i+1>Nx and ijkMesh_.symXOut())
	       il=iP+1-i;
	     else if (iP+i<0 and ijkMesh_.symXIn())
	       il=iP-1-i;
	     else
	       il=(iP+i+Nx)%Nx;
	     for (int j=-nMax[1];j<nMax[1]+1;j++)
	       {
		 if (jP+j+1>Ny and ijkMesh_.symYOut()) 
		   jl=jP+1-j;
		 else if (jP+j<0 and ijkMesh_.symYIn())
		   jl=jP-1-j;
		 else
		   jl=(jP+j+Ny)%Ny;
		 label ijkG=ijkMesh_.ijk1(il,jl,kl); 
		 label gblIdx=globalIds[ijkG];
		 if (dir==0)
		   H[j+1]+=phiIJK[gblIdx];
		 else
		   H[i+1]+=phiIJK[gblIdx];
	       }
	   }
       }
   }
  else
    {
      H.setSize(9,0.0);
      for (int i=-nMax[0];i<nMax[0]+1;i++)
	{
	  if (iP+i+1>Nx and ijkMesh_.symXOut())
	    il=iP+1-i;
	  else if (iP+i<0 and ijkMesh_.symXIn())
	    il=iP-1-i;
	  else
	    il=(iP+i+Nx)%Nx;
	  for (int j=-nMax[1];j<nMax[1]+1;j++)
	    {
	      if (jP+j+1>Ny and ijkMesh_.symYOut()) 
		jl=jP+1-j;
	      else if (jP+j<0 and ijkMesh_.symYIn())
		jl=jP-1-j;
	      else
		jl=(jP+j+Ny)%Ny;
	      for (int k=-nMax[2];k<nMax[2]+1;k++)
		{
		  if (kP+k+1>Nz and ijkMesh_.symZOut()) 
		    kl=kP+1-k;
		  else if (kP-k<0 and ijkMesh_.symZIn())
		    kl=kP-1-k;
		  else
		    kl=(kP+k+Nz)%Nz;
		  label ijkG=ijkMesh_.ijk1(il,jl,kl); 
		  label gblIdx=globalIds[ijkG];
		  if (dir==0)
		    H[j+1+3*(k+1)]+=phiIJK[gblIdx];
		  else if (dir==1)
		    H[i+1+3*(k+1)]+=phiIJK[gblIdx];
		  else
		    H[i+1+3*(j+1)]+=phiIJK[gblIdx];
		}
	    }
	}
    }

  return -H*dl_;
  
}

void Foam::heightFunction::calculateK()
{
  Info<<"Calculating the curvature with the height function method..."<<endl;
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

  // no solid boundaries yet
  //const volVectorField gradAlpha(fvc::grad(alpha1_, "nHat"));
  //surfaceVectorField gradAlphaf(fvc::interpolate(gradAlpha));
  //surfaceVectorField nHatfv(gradAlphaf/(mag(gradAlphaf) + deltaN_));
  // //const volVectorField nHat(gradAlpha/(mag(gradAlpha) + deltaN_));
  //correctContactAngle(nHatfv.boundaryFieldRef(), gradAlphaf.boundaryFieldRef());

  forAll(K_,celli)
    K_[celli]=0;


  reconstructionSchemes& surf =
    mesh.lookupObjectRef<reconstructionSchemes>("reconstructionScheme");

  surf.reconstruct(false);
  //surf.reconstruct(true);

  const boolList& interfaceCells = surf.interfaceCell();
    
  const volVectorField& faceCentre = surf.centre();
  const volVectorField& faceNormal = surf.normal();

  Vector<label> stencilSize(nMax_,nMax_,nMax_);
  ijkMesh_.getZoneField(interfaceCells,alphaIJK,alpha1_,stencilSize);
  
  forAll(interfaceCells,celli)
    {
      vector n = faceNormal[celli];
      if(interfaceCells[celli] && mag(n)>0)
	//vector n = nHat[celli];
      //RDF_.correctBoundaryConditions();
      //RDF_.markCellsNearSurf(interfaceCells,1);
      //const boolList& nextToInterface =RDF_.nextToInterface(); 
      //if(interfaceCells[celli] or nextToInterface[celli])
	{
	  label dir=-1;
	  if (mag(n.x())>0.0 and !ijkMesh_.isEmpty().x())
	    dir=0;
	  if (mag(n.y())>mag(n.x()) and !ijkMesh_.isEmpty().y())
	    dir=1;
	  if (mag(n.z())>mag(n.y()) and mag(n.z())>mag(n.x()) and !ijkMesh_.isEmpty().z())
	    dir=2;
	  if (dir==-1)
	    {
	      Info<<"n="<<n<<endl;
	      FatalErrorInFunction
	      << "Normal direction cannot be identified."
	     << abort(FatalError);
	    }

	  if (boundaryCells_[celli])
	    {
	      //Info<<"C[celli]="<<C[celli]<<endl;
	      printf("Cell at the boundary: proc=%d, x=%f, y=%f, z=%f\n",Pstream::myProcNo(),C[celli].x(),C[celli].y(),C[celli].z());
	      /*FatalErrorInFunction
		<< "Interface is at a non-cyclic cellSet or domain boundary. Non-cyclic boundaries is not supported at the moment."
		<< abort(FatalError);*/
	      continue;
	    }

	  //stencil_.setStencil(alphaIJK,ijkG);
	  //K_[celli]=stencil_.calcCurvature(dir);

	  Vector<label> ijkG=ijkMesh_.ijk3(celli);
	  List<scalar> H=calculateHeights(alphaIJK,ijkG,dir); 

	  if (is2D_)
	    {
	      scalar Ht=(H[2]-H[0])/2.0/dl_;
	      scalar Htt=(H[2]-2.0*H[1]+H[0])/dl_/dl_;
	      K_[celli]=Htt/Foam::pow(1.0+Ht*Ht,1.5);
	    }
	  else
	    {
	      double H2[3][3]={{0.0,0.0,0.0},{0.0,0.0,0.0},{0.0,0.0,0.0}}; 
	      for (int i=-1;i<2;i++)//normal 
		{
		  for (int j=-1;j<2;j++)//tangent
		    {
		      H2[i+1][j+1]=H[i+1+3*(j+1)]; 
		    }
		}
	      scalar Ht=(H2[2][1]-H2[0][1])/2.0/dl_;
	      scalar Hb=(H2[1][2]-H2[1][0])/2.0/dl_;
	      scalar Htt=(H2[2][1]-2.0*H2[1][1]+H2[0][1])/dl_/dl_;
	      scalar Hbb=(H2[1][2]-2.0*H2[1][1]+H2[1][0])/dl_/dl_;
	      scalar Htb=(H2[2][2]-H2[2][0]-H2[0][2]+H2[0][0])/4./dl_/dl_;
	      K_[celli]=(Htt+Hbb+Htt*Hb*Hb+Hbb*Ht*Ht-2.0*Htb*Ht*Hb)/Foam::pow(1.0+Ht*Ht+Hb*Hb,1.5);      
	    }	  	      
	}      
    }

#include "../interpolateCutCellsToFaces.H"
}



// ************************************************************************* //
