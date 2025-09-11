forAll(interfaceCells,iCell)
{
  if(interfaceCells[iCell])
    {
      //neiInterface[iCell]=true;
      if (alpha1_[iCell]<interfaceTol_ or alpha1_[iCell]>1-interfaceTol_) //mag(faceCentre)>
	{
	  K_[iCell]=0.0;
	}
    }
}

K_.correctBoundaryConditions();
Kf_=fvc::interpolate(K_);

//const surfaceVectorField& Cf = mesh.Cf();
const Foam::labelList& nei=alpha1_.mesh().faceNeighbour();
const Foam::labelList& own=alpha1_.mesh().faceOwner();
// If a face is between a cutcell and empty cell, its curvature value is that of cutcell.
forAll(Kf_,iFace)
{
  label iCellO=own[iFace];
  label iCellN=nei[iFace];
  scalar K1=K_[own[iFace]];      
  scalar K2=K_[nei[iFace]];
  if (Kf_[iFace]!=0 and K1*K2==0.0)
{
    if(Foam::mag(K1)>Foam::mag(K2))
	Kf_[iFace]=K1;
    else
	Kf_[iFace]=K2;      
  }
}

//Interior faces: if exactly one side is a cut cell, use the cut-cell curvature
//const label nIntFaces = mesh.nInternalFaces();
//const labelList& own = mesh.faceOwner();
//const labelList& nei = mesh.faceNeighbour();

//for (label f = 0; f < nIntFaces; ++f)
//{
//  const label cO = own[f];
//  const label cN = nei[f];

//  const scalar aO = alpha1_[cO];
//  const scalar aN = alpha1_[cN];

//  const bool cutO = (aO > interfaceTol_ && aO < 1.0 - interfaceTol_);
//  const bool cutN = (aN > interfaceTol_ && aN < 1.0 - interfaceTol_);

//  if (cutO ^ cutN)   // exactly one cut
//  {
//      Kf_[f] = cutO ? K_[cO] : K_[cN];
//  }
    // else: both cut -> keep linear Kf_[f]; both pure -> Kf_[f] should already be ~0
//}

// If a boundary face has non-zero curvature, we assume it is on a periodic boundary -- must generalize later!
// We re-apply linear interpolation if neigbouring cells of the periodic face are interface cells. If not, only cutcell-value is used. 

scalar thr=1e-14/dl_;
forAll(mesh.boundaryMesh(),iPatch)
{
  word pType=Kf_.boundaryField()[iPatch].type();
  if(pType!="empty")
    {	  
      const polyPatch& cPatch = mesh.boundaryMesh()[iPatch];
      //label iFaceStart = cPatch.start();
      const labelUList& faceCells = cPatch.faceCells();
      forAll(faceCells,iFace)
	{
	  if (Kf_.boundaryFieldRef()[iPatch][iFace]!=0)
	    {
	      //Kf_.boundaryFieldRef()[iPatch][iFace]=0;
	      //const point Cfp=Cf.boundaryField()[iPatch][iFace];

	      label iCell=faceCells[iFace];
	      scalar K=K_[iCell];	    
	      //if Kf=K/2. or K=0, it means one neighbouring cell is empty, and interpolation was halved (K/2). Correction by doubling:
	      if (Foam::mag(2.*Kf_.boundaryFieldRef()[iPatch][iFace]-K)<thr or K<thr)
		//Kf_.boundaryFieldRef()[iPatch][iFace]=K;
		Kf_.boundaryFieldRef()[iPatch][iFace]*=2.;

	      //Info<<"K_["<<celli<<"]="<<K_[celli]<<endl;
	    }
	}
    }
}  

 

