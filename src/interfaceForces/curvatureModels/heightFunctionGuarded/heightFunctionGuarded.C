/*---------------------------------------------------------------------------*\
            Copyright (c) 2026, Asim Onder
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

#include "heightFunctionGuarded.H"
#include "addToRunTimeSelectionTable.H"

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
#include "DynamicList.H"
#include "reconstructedDistanceFunction.H"
#include "rdfCurvature.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(heightFunctionGuarded, 0);
    addToRunTimeSelectionTable(curvatureModel,heightFunctionGuarded, components);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::heightFunctionGuarded::heightFunctionGuarded
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
    ijkMesh_(ijkZone::New(alpha1.mesh())),
    globalNumbering_(ijkMesh_.globalNumbering()),
    nMax_
    (
        dict.subOrEmptyDict("heightFunctionParams")
            .lookupOrDefault<label>("nMax",3)
    ),
    nMaxExtend_
    (
        dict.subOrEmptyDict("heightFunctionParams")
            .lookupOrDefault<label>("nMaxExtend",8)
    ),
    tolPure_
    (
        dict.subOrEmptyDict("heightFunctionParams")
            .lookupOrDefault<scalar>("tolPure",1e-9)
    ),
    interfaceTol_(1e-12),
    boundaryCells_(mesh_.nCells(),false)
{

  Info<<"Curvature model: heightFunctionGuarded."
      <<" nMax="<<nMax_<<" nMaxExtend="<<nMaxExtend_
      <<" tolPure="<<tolPure_<<endl;

  if (nMaxExtend_ < nMax_)
  {
      nMaxExtend_ = nMax_;
  }

  // margin matches the maximum window reach
  ijkMesh_.markBoundaryCells(boundaryCells_,nMaxExtend_);

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
}

// * * * * * * * * * * * * * * Public Access Member Functions  * * * * * * * //

void Foam::heightFunctionGuarded::correctContactAngle
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


// index rules identical to heightFunction::calculateHeights, including
// the z-axis symmetry-in condition as it stands there

Foam::label Foam::heightFunctionGuarded::wrapX(label b, label o) const
{
    const label Nx=ijkMesh_.Nx();
    if (b+o+1>Nx and ijkMesh_.symXOut()) return b+1-o;
    if (b+o<0 and ijkMesh_.symXIn())     return b-1-o;
    return (b+o+Nx)%Nx;
}

Foam::label Foam::heightFunctionGuarded::wrapY(label b, label o) const
{
    const label Ny=ijkMesh_.Ny();
    if (b+o+1>Ny and ijkMesh_.symYOut()) return b+1-o;
    if (b+o<0 and ijkMesh_.symYIn())     return b-1-o;
    return (b+o+Ny)%Ny;
}

Foam::label Foam::heightFunctionGuarded::wrapZ(label b, label o) const
{
    const label Nz=ijkMesh_.Nz();
    if (b+o+1>Nz and ijkMesh_.symZOut()) return b+1-o;
    if (b-o<0 and ijkMesh_.symZIn())     return b-1-o;
    return (b+o+Nz)%Nz;
}


Foam::heightFunctionGuarded::tryStatus
Foam::heightFunctionGuarded::tryHeights
(
    const Map<scalar>& phiIJK,
    const Vector<label>& ijk,
    const label dir,
    const label m,
    List<scalar>& H
) const
{
  const Vector<bool>& isEmpty=ijkMesh_.isEmpty();

  const label iP=ijk.x();
  const label jP=ijk.y();
  const label kP=ijk.z();

  // cross axes of dir: 3D uses both, 2D uses the one in-plane axis
  label axA=-1, axB=-1;
  if (dir==0)      { axA=1; axB=2; }
  else if (dir==1) { axA=0; axB=2; }
  else             { axA=0; axB=1; }
  if (is2D_)
  {
      if (isEmpty.x())      { axA = (dir==1 ? 2 : 1); axB=-1; }
      else if (isEmpty.y()) { axA = (dir==0 ? 2 : 0); axB=-1; }
      else                  { axA = (dir==0 ? 1 : 0); axB=-1; }
  }

  const label nSlots = is2D_ ? 3 : 9;
  H.setSize(nSlots);
  H = 0.0;

  label orient=0;

  const label pLo=-1, pHi=1;
  const label qLo = is2D_ ? 0 : -1;
  const label qHi = is2D_ ? 0 : 1;

  for (label p=pLo; p<=pHi; ++p)
  {
      for (label q=qLo; q<=qHi; ++q)
      {
          scalar sum=0.0, aPrev=0.0, aLo=0.0, aHi=0.0;
          label nCross=0;

          for (label w=-m; w<=m; ++w)
          {
              label di=0, dj=0, dk=0;

              if (dir==0) di=w; else if (dir==1) dj=w; else dk=w;
              if (axA==0) di+=p; else if (axA==1) dj+=p; else if (axA==2) dk+=p;
              if (axB==0) di+=q; else if (axB==1) dj+=q; else if (axB==2) dk+=q;

              const label il = isEmpty.x() ? 0 : wrapX(iP,di);
              const label jl = isEmpty.y() ? 0 : wrapY(jP,dj);
              const label kl = isEmpty.z() ? 0 : wrapZ(kP,dk);

              const label gblIdx=ijkMesh_.lookupGlobalId(il,jl,kl);
              if (!phiIJK.found(gblIdx)) return MISSING;

              const scalar a=phiIJK[gblIdx];
              if (w==-m) aLo=a;
              if (w== m) aHi=a;
              if (w>-m and (aPrev-0.5)*(a-0.5)<0.0) ++nCross;
              aPrev=a;
              sum+=a;
          }

          // guard 2: a single interface in the column
          if (nCross>1) return IMPURE;

          // guard 1: pure ends, one shared orientation
          label o=0;
          if (aLo>1.0-tolPure_ and aHi<tolPure_)      o= 1;
          else if (aLo<tolPure_ and aHi>1.0-tolPure_) o=-1;
          else return IMPURE;
          if (orient==0) orient=o;
          else if (orient!=o) return IMPURE;

          const label slot = is2D_ ? p+1 : (p+1) + 3*(q+1);
          H[slot] = -sum*dl_;
      }
  }

  return OK;
}


Foam::scalar Foam::heightFunctionGuarded::kappaFromHeights
(
    const List<scalar>& H
) const
{
  if (is2D_)
  {
      const scalar Ht=(H[2]-H[0])/2.0/dl_;
      const scalar Htt=(H[2]-2.0*H[1]+H[0])/dl_/dl_;
      return Htt/Foam::pow(1.0+Ht*Ht,1.5);
  }

  scalar H2[3][3];
  for (label i=-1;i<2;i++)
  {
      for (label j=-1;j<2;j++)
      {
          H2[i+1][j+1]=H[i+1+3*(j+1)];
      }
  }

  const scalar Ht=(H2[2][1]-H2[0][1])/2.0/dl_;
  const scalar Hb=(H2[1][2]-H2[1][0])/2.0/dl_;
  const scalar Htt=(H2[2][1]-2.0*H2[1][1]+H2[0][1])/dl_/dl_;
  const scalar Hbb=(H2[1][2]-2.0*H2[1][1]+H2[1][0])/dl_/dl_;
  const scalar Htb=(H2[2][2]-H2[2][0]-H2[0][2]+H2[0][0])/4./dl_/dl_;
  return
      (Htt+Hbb+Htt*Hb*Hb+Hbb*Ht*Ht-2.0*Htb*Ht*Hb)
     /Foam::pow(1.0+Ht*Ht+Hb*Hb,1.5);
}


void Foam::heightFunctionGuarded::calculateK()
{
  const fvMesh& mesh = alpha1_.mesh();
  const volVectorField& C = mesh.C();
  Map<scalar> alphaIJK;

  forAll(K_,celli)
    K_[celli]=0;

  reconstructionSchemes& surf =
    mesh.lookupObjectRef<reconstructionSchemes>("reconstructionScheme");

  const boolList& interfaceCells = surf.interfaceCell();
  const volVectorField& faceNormal = surf.normal();

  Vector<label> stencilSize(nMaxExtend_,nMaxExtend_,nMaxExtend_);
  ijkMesh_.getZoneField(interfaceCells,alphaIJK,alpha1_,stencilSize);

  const Vector<bool>& isEmpty=ijkMesh_.isEmpty();

  boolList kValid(mesh_.nCells(),false);
  DynamicList<label> needNbr;
  DynamicList<label> fallbackCells;

  label nHFCells=0, nStd=0, nDirRetry=0, nGrown=0, nBoundaryFallback=0;

  forAll(interfaceCells,celli)
    {
      const vector n = faceNormal[celli];
      if(interfaceCells[celli] && mag(n)>0)
        {
          nHFCells++;

          if (boundaryCells_[celli])
            {
              fallbackCells.append(celli);
              nBoundaryFallback++;
              continue;
            }

          // directions in decreasing order of the normal component,
          // empty axes excluded
          scalar nc[3]={mag(n.x()),mag(n.y()),mag(n.z())};
          if (isEmpty.x()) nc[0]=-1;
          if (isEmpty.y()) nc[1]=-1;
          if (isEmpty.z()) nc[2]=-1;

          DynamicList<label> dirs;
          for (label d=0; d<3; ++d)
              if (nc[d]>=0) dirs.append(d);
          // selection sort, descending, stable (x before y before z)
          for (label a=0; a<dirs.size(); ++a)
          {
              label best=a;
              for (label b=a+1; b<dirs.size(); ++b)
                  if (nc[dirs[b]]>nc[dirs[best]]) best=b;
              if (best!=a) Swap(dirs[a],dirs[best]);
          }

          const Vector<label> ijkG=ijkMesh_.ijk3(celli);
          bool found=false;

          forAll(dirs,dd)
            {
              const label dir=dirs[dd];
              for (label m=nMax_; m<=nMaxExtend_; ++m)
                {
                  List<scalar> H;
                  const tryStatus st=tryHeights(alphaIJK,ijkG,dir,m,H);
                  if (st==OK)
                    {
                      K_[celli]=kappaFromHeights(H);
                      kValid[celli]=true;
                      found=true;
                      if (dd==0 and m==nMax_) nStd++;
                      else if (m==nMax_)      nDirRetry++;
                      else                    nGrown++;
                      break;
                    }
                  if (st==MISSING) break;   // larger windows miss too
                }
              if (found) break;
            }

          if (!found) needNbr.append(celli);
        }
    }

  // --- neighbour-average pass: valid pass-1 curvature only ---------------
  label nNbr=0;
  if (returnReduce(needNbr.size(), sumOp<label>()) > 0)
  {
      volScalarField Ktmp
      (
          IOobject
          (
              "KtmpGuarded",
              mesh_.time().timeName(),
              mesh_,
              IOobject::NO_READ,
              IOobject::NO_WRITE
          ),
          K_
      );
      forAll(Ktmp,celli)
          if (!kValid[celli]) Ktmp[celli]=0.0;

      Map<scalar> kIJK;
      Vector<label> one(1,1,1);
      ijkMesh_.getZoneField(interfaceCells,kIJK,Ktmp,one);

      forAll(needNbr,ii)
        {
          const label celli=needNbr[ii];
          const Vector<label> ijkG=ijkMesh_.ijk3(celli);

          scalar sk=0.0; label a=0;
          for (label di=-1; di<2; ++di)
          for (label dj=-1; dj<2; ++dj)
          for (label dk=-1; dk<2; ++dk)
            {
              if (di==0 and dj==0 and dk==0) continue;
              if (isEmpty.x() and di!=0) continue;
              if (isEmpty.y() and dj!=0) continue;
              if (isEmpty.z() and dk!=0) continue;

              const label il = isEmpty.x() ? 0 : wrapX(ijkG.x(),di);
              const label jl = isEmpty.y() ? 0 : wrapY(ijkG.y(),dj);
              const label kl = isEmpty.z() ? 0 : wrapZ(ijkG.z(),dk);

              const label gblIdx=ijkMesh_.lookupGlobalId(il,jl,kl);
              if (!kIJK.found(gblIdx)) continue;
              const scalar v=kIJK[gblIdx];
              if (mag(v)>VSMALL) { sk+=v; ++a; }
            }
          if (a>0)
            {
              K_[celli]=sk/a;
              kValid[celli]=true;
              nNbr++;
            }
          else
            {
              fallbackCells.append(celli);
            }
        }
  }

  // --- rdfCurvature, the last resort (single-path octree) ----------------
  label nRDF=fallbackCells.size();
  label nFallbackGlobal = nRDF;
  reduce(nFallbackGlobal, sumOp<label>());

  if (nFallbackGlobal > 0)
    {
      tmp<volScalarField> tKRDF = rdfCurvature::computeKRDF(mesh, alpha1_);
      const volScalarField& KRDF = tKRDF();

      forAll(fallbackCells, i)
        {
          K_[fallbackCells[i]] = KRDF[fallbackCells[i]];
        }
    }

  reduce(nHFCells, sumOp<label>());
  reduce(nStd, sumOp<label>());
  reduce(nDirRetry, sumOp<label>());
  reduce(nGrown, sumOp<label>());
  reduce(nNbr, sumOp<label>());
  reduce(nRDF, sumOp<label>());
  reduce(nBoundaryFallback, sumOp<label>());

  Info<< "heightFunctionGuarded: " << nHFCells << " interface cells: "
      << nStd << " standard, " << nDirRetry << " direction-retry, "
      << nGrown << " grown-window, " << nNbr << " neighbour-average, "
      << nRDF << " rdf (" << nBoundaryFallback << " zone-boundary)" << endl;

  if (is2D_)
  {
    #include "../interpolateCutCellsToFaces2D.H"
  }
  else
  {
    #include "../interpolateCutCellsToFaces3D.H"
  }
}

// ************************************************************************* //
