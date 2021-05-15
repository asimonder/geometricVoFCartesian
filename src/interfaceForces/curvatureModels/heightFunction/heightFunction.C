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
    deltaN_
    (
        "deltaN",
        1e-8/pow(average(alpha1.mesh().V()), 1.0/3.0)
     ),
    mesh_(alpha1.mesh()),
    RDF_(alpha1.mesh()),
    ijkMesh_(alpha1.mesh()),
    globalNumbering_(ijkMesh_.globalNumbering()),
    nMax_(3),
    boundaryCells_(mesh_.nCells(),false)
{
  
  ijkMesh_.markBoundaryCells(boundaryCells_,nMax_);
  Info<<"Done setting up ijkMesh for the height function method..."<<endl;
  
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


void Foam::heightFunction::calculateK()
{
  Info<<"Calculting the curvature with the heigh function method..."<<endl;
  const label iMax=nMax_;
  label jMax,kMax=0;
  //const scalar nEls=2.0*iMax+1;
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

  forAll (K_,celli)
    K_[celli]=0;

  // where to put this one?
  const volVectorField gradAlpha(fvc::grad(alpha1_, "nHat"));
  surfaceVectorField gradAlphaf(fvc::interpolate(gradAlpha));
  surfaceVectorField nHatfv(gradAlphaf/(mag(gradAlphaf) + deltaN_));
  correctContactAngle(nHatfv.boundaryFieldRef(), gradAlphaf.boundaryFieldRef());
    
  reconstructionSchemes& surf =
    mesh.lookupObjectRef<reconstructionSchemes>("reconstructionScheme");

  surf.reconstruct(false);

  const boolList& interfaceCells = surf.interfaceCell();

    
  const volVectorField& faceCentre = surf.centre();
  const volVectorField& faceNormal = surf.normal();

  Vector<label> stencilSize(iMax,iMax,iMax);
  ijkMesh_.getZoneField(interfaceCells,alphaIJK,alpha1_,stencilSize);

  scalar dN=dx;
  scalar dT=dz;
  if (Nz==1)
    dT=dy;

  //Info<<"2: Calculting the curvature with the heigh function method..."<<endl;
  forAll(interfaceCells,celli)
    {
      vector n = faceNormal[celli];
      if(interfaceCells[celli] && mag(n)>0)
	{
	  label dir=0;
	  scalar nHat=0.0;
	  if (mag(n.x())>0.0)
	    nHat=n.x()/mag(n.x());
	  const point cc = C[celli];
	  
	  if (mag(n.y())>mag(n.x()) && Ny>1)
	    {
	      dir=1;
	      nHat=n.y()/mag(n.y());
	      dN=dy;dT=dx;
	    }
	  if (mag(n.z())>mag(n.y()) && mag(n.z())>mag(n.x()) && Nz>1)
	    {
	      dir=2;
	      nHat=n.z()/mag(n.z());
	      dN=dz;dT=dx;
	    }

	  if (boundaryCells_[celli])
	    {
	      FatalErrorInFunction
		<< "Interface is at a non-cyclic cellSet or domain boundary. Non-cyclic boundaries is not supported at the moment."
		<< abort(FatalError);
	    }
	      
	  label i=round((cc.x()-Pmin.x())/dx);
	  label j=round((cc.y()-Pmin.y())/dy);
	  label k=round((cc.z()-Pmin.z())/dz);
	  if (Ny==1)
	    j=0;
	  if (Nz==1)
	    k=0;

	  double H[3]={0.0,0.0,0.0}; 
	  label il,jl,kl;
	  for (int iN=-iMax;iN<iMax+1;iN++)
	    {
	      for (int iT=-1;iT<2;iT++)
		{
		  if (dir==0)
		    {
		      il=(i+iN)%Nx;
		      if (Ny==1)
			{
			  kl=(k+iT);
			  jl=0;
			}
		      else
			{
			  jl=(j+iT);
			  kl=0;
			}
		    }
		  else if (dir==1)
		    {
		      il=(i+iT)%Nx;jl=(j+iN);kl=0;
		    }
		  else
		    {
		      il=(i+iT)%Nx;kl=(k+iN);jl=0;
		    }			    			   
		  label ijk=il+Nx*jl+Nx*Ny*kl;	  
		  label gblIdx=globalIds[ijk];
		  H[iT+1]+=alphaIJK[gblIdx]*dN;
		}
	    }
	  scalar Ht=(H[2]-H[0])/2.0/dT;
	  scalar Htt=(H[2]-2.0*H[1]+H[0])/dT/dT;
	  scalar kappa=(Htt)/Foam::pow(1.0+Ht*Ht,1.5);
	  K_[celli]=-kappa; 
	}
      
    }

  //Info<<"3: Calculting the curvature with the heigh function method..."<<endl;  
  label neiMax=1;
  //this one can cause problems at cyclic BC!
  RDF_.markCellsNearSurf(interfaceCells,neiMax);
  const boolList& nextToInterface =RDF_.nextToInterface(); 
  //boolList atBoundary=ijkMesh_.markBoundaryCells(1);
  
  boolList neiInterface(mesh.nCells(),false);

  forAll(neiInterface,celli)
    {
      if (nextToInterface[celli] and !interfaceCells[celli]) // and !boundaryCells_[celli])
	  neiInterface[celli]=true;
    }

  stencilSize=Vector<label>(neiMax+1,neiMax+1,neiMax+1);
  ijkMesh_.getZoneField(neiInterface,faceCentreIJK,faceCentre,stencilSize);
  ijkMesh_.getZoneField(neiInterface,curvIJK,K_,stencilSize);

  //Info<<"4: Calculting the curvature with the heigh function method..."<<endl;
  
  neiMax+=1;
  forAll(neiInterface,celli)
    {
      if (neiInterface[celli])
	{
	  const point cc = C[celli];

	  label i=round((cc.x()-Pmin.x())/dx);
	  label j=round((cc.y()-Pmin.y())/dy);
	  label k=round((cc.z()-Pmin.z())/dz);
	  if (Ny==1)
	    {
	      j=0;
	      jMax=0;
	      kMax=neiMax;
	    }
	  
	  if (Nz==1)
	    {
	      k=0;
	      kMax=0;
	      jMax=neiMax;
	    }
	  label gblIdMin=-1;
	  scalar distMin=GREAT;
	  for (int ii=-neiMax;ii<neiMax+1;ii++)
	    {
	      for (int jj=-jMax;jj<jMax+1;jj++)
		{
		  for (int kk=-kMax;kk<kMax+1;kk++)
		    {
		      label il=(i+ii)%Nx;
		      label kl=(k+kk); 
		      label jl=(j+jj);
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

  // Info<<"5: Calculting the curvature with the heigh function method..."<<endl;
}



// ************************************************************************* //
