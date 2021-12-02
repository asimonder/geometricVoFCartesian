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

Application
    generateCurvatureData

Description


Author
   Asim Onder

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

#include "OFstream.H"
#include "autoPtr.H"
#include "reconstructionSchemes.H"
#include "implicitFunction.H"
#include "cutCellImpFunc.H"
#include "cutCellIso.H"
#include "reconstructionError.H"
#include "writeFile.H"
//#include "isoCutCell.H"

#include "interfaceForces.H"
//#include "OBJstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


void setAlpha
(
    const fvMesh& mesh,
    const dictionary& setAlphaFieldDict,
    volScalarField& alpha1
)
{

  Foam::autoPtr<implicitFunction> func(Foam::implicitFunction::New
				       (
					setAlphaFieldDict.get<word>("type"),
					setAlphaFieldDict
					));	    
  word functionType (setAlphaFieldDict.lookup("type"));
  bool is2D=(functionType == "cylinder");
  const vector origin=setAlphaFieldDict.get<vector>("origin");
  scalar radius=setAlphaFieldDict.get<scalar>("radius");
  const cellList& cells=mesh.cells();
  const faceList& faces=mesh.faces();
  const  surfaceVectorField&  Cf = mesh.Cf();
  scalar dx=0.0;
  label Nx=200;
  label Ny=Nx;
  if (is2D)
    Ny=1;
  label Nz=Nx;

  //only works for uniform hexa cells aligned with Cartesian directions
  forAll(alpha1,iCell)
    {
      vector Pmin=vector::uniform(VGREAT);
      vector Pmax=vector::uniform(-VGREAT);
      bool isWater=true;
      bool isAir=true;
      for (int iF=0;iF<6;iF++)
	{
	  label f1=cells[iCell][iF];
	  Pmin.x()=Foam::min(Cf[f1].x(),Pmin.x());
	  Pmin.y()=Foam::min(Cf[f1].y(),Pmin.y());
	  Pmin.z()=Foam::min(Cf[f1].z(),Pmin.z());
	  Pmax.x()=Foam::max(Cf[f1].x(),Pmax.x());
	  Pmax.y()=Foam::max(Cf[f1].y(),Pmax.y());
	  Pmax.z()=Foam::max(Cf[f1].z(),Pmax.z());
	}

      dx=(Pmax.x()-Pmin.x());
      for (int i=0;i<2;i++)
	{
	  for (int j=0;j<2;j++)
	    {
	      for (int k=0;k<2;k++)
		{
		  vector Pc=Pmin+vector(dx*i,dx*j,dx*k);
		  if (mag(Pc-origin)<radius)
		      isAir=false;
		  if (mag(Pc-origin)>radius)
		      isWater=false;
		}
	    }
	}
     
      //Info<<"dx="<<dx<<endl;
      if (isWater)
	alpha1[iCell]=1.0;
      else if (isAir)
	alpha1[iCell]=0.0;
      else
	{
	  dx=(Pmax.x()-Pmin.x())/Nx;
	  vector P0=Pmin+vector(dx/2.,dx/2.,dx/2.);
	  scalar sum=0;
	  for (int i=0;i<Nx;i++)
	    {
	      for (int j=0;j<Ny;j++)
		{
		  for (int k=0;k<Nz;k++)
		    {
		      vector Pc=P0+vector(dx*i,dx*j,dx*k);
		      if (mag(Pc-origin)>radius)
			sum+=1.0;
		    }
		}
	    }
	  alpha1[iCell]=sum/(Nx*Ny*Nz);
	}
    }

  alpha1.correctBoundaryConditions();
}





int main(int argc, char *argv[])
{
//    #include "addRegionOption.H"
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createFields.H"
  
    ijkZone ijkMesh(mesh);
    const label NCells=ijkMesh.Nx();
    //vector Pmin=ijkMesh.Pmin();  
    //vector Pmax=ijkMesh.Pmax();  


    #include "createTimeControls.H"
    #include "readTimeControls.H"
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"

    writeFile writeObj(mesh, "curvatureData", "curvDataSet");
      
    Random rndCentre(1234567);

    autoPtr<reconstructionSchemes> surf =
        reconstructionSchemes::New(alpha1,phi,U,fvSolutionDict);
    interfaceForces surfForces(alpha1,U,transportProperties);


    label nIter = fvSolutionDict.get<label>("nIter");
    word setAlphaMethod = fvSolutionDict.get<word>("setAlphaMethod");
    if (setAlphaMethod != "cutCellImpFunc" && setAlphaMethod != "cutCellIso")
    {
        FatalError  << "valid choice are only cutCellImpFunc or cutCellIso"
                    << abort(FatalError);
    }

    word functionType (setAlphaFieldDict.lookup("type"));
    bool twoDim = (functionType == "cylinder");
    Info << "twoDim = " << twoDim << endl;
    Info << "functionType = " << functionType << endl;

    scalar recTime = 0;
    //vector centreMin = Pmin+0.2*Pmin;
    //vector centreMax = Pmax+0.2*Pmax;;
    //scalar Rmin=0.001;
    //scalar dR=0.001;
    //label nIter=5;
    label nMax=1;
    boolList interfaceCells(mesh.nCells(),false);
    uniformStencil stencil(mesh,ijkMesh,nMax);
    //scalar thr=0.000001;
    const scalar interfaceTol=setAlphaFieldDict.get<scalar>("interfaceTol");

    
    FILE *writeTraining;
    FILE *writeNormals;
    FILE *writeCurvatures;
    string strNc=std::to_string(NCells);
    string dictPath=runTime.path()/"trainingData_N"+strNc+".dat";
    writeTraining = fopen (dictPath.c_str(),"w");
    dictPath=runTime.path()/"normalData__N"+strNc+"_"+reconScheme+".dat";
    writeNormals = fopen (dictPath.c_str(),"w");
    dictPath=runTime.path()/"curvatureData_N"+strNc+"_"+curvModel+"_"+reconScheme+".dat";
    writeCurvatures = fopen (dictPath.c_str(),"w");
    

    
    while (runTime.run())
    {

      Foam<implicitFunction> func(Foam::implicitFunction::New
				  (
				   setAlphaFieldDict.get<word>("type"),
				   setAlphaFieldDict
				   ));

      scalarField f(mesh.nPoints(),0.0);
      forAll(f,pI)
	{
	  f[pI]=func->value(mesh.points()[pI]);
	}
      
      cutCellImpFunc cutCell(mesh,f,func.ref());
      setAlpha(cutCell,mesh,setAlphaFieldDict,alpha1);
	      
      surf->reconstruct();	      


      volScalarField curv = mesh.lookupObjectRef<volScalarField>("K_");

	      //select interfacecells
      interfaceCells=false;
      if (alpha1,iCell)
	{
	  if (alpha[iCell]>thr and alpha[iCell]<1-thr)
	    interfaceCells[iCell]=true;
	}
      
      Map<scalar> phiIJK;
      Vector<label> stencilSize(nMax,nMax,nMax);
      ijkMesh_.getZoneField(interfaceCellw,phiIJK,phi,stencilSize);
      
      List<scalar> A;
      forAll(interfaceCells, i)		
	    {
	      if(interfaceCells[iCell])
		{
		  stencil_.setStencil(phiIJK,ijkMesh_.ijk3(iCell));
		  A=stencil_getStencil();
		}
	    }
		  
	  //write to a file
	  if (Pstream::master())
	    {
	      //writeTime(file());
	      
	      writeObj.file()<< curv;
	      
	      forAll(A,i)
		writeObj.file()<< A[i];
	      
	      writeObj.file()<< endl;
	    }
	  
	}

      runTime.write();

    }


    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
