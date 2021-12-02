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
    generateJONSWAPData

Description


Author
   Asim Onder

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

#include "OFstream.H"
#include "autoPtr.H"
#include "reconstructionSchemes.H"
#include "waveFunction.H"
#include "writeFile.H"
//#include "isoCutCell.H"

#include "interfaceForces.H"
#include "ijkZone.H"
#include "uniformStencil.H"
//#include "randomWaveFieldImplicitFunction.H"
//#include "OBJstream.H"
#include <ctime>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void setAlpha
(
    const fvMesh& mesh,
    //const dictionary& setAlphaFieldDict,
    Foam::autoPtr<waveFunction>& func,
    volScalarField& alpha1
)
{

  const cellList& cells=mesh.cells();
  const faceList& faces=mesh.faces();
  const  surfaceVectorField&  Cf = mesh.Cf();
  label Nx=1;
  label Ny=1;
  label iC=0;
  forAll(alpha1,iCell)
    {
      iC+=1;
      vector Pmin=vector::uniform(VGREAT);
      vector Pmax=vector::uniform(-VGREAT);
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

      scalar dxC=(Pmax.x()-Pmin.x());

      scalar dx=dxC/Nx;

      //printf("proc=%d: iCell=%d/%d, dx=%f, dxC=%f\n",Pstream::myProcNo(),iC, mesh.nCells(),dx,dxC);
      
      scalar V=dxC*dxC*dxC;
      vector P0=Pmin+vector(dx/2.,dx/2.,0);
      scalar sum=0;
      scalar Hmax=-VGREAT;
      scalar Hmin=VGREAT;
      for (int i=0;i<Nx;i++)
	    {
	      for (int j=0;j<Ny;j++)
		{
		  vector Pc=P0+vector(dx*i,dx*j,0);
		  scalar val=func->value(Pc);
		  Hmax=Foam::max(Hmax,val);
		  Hmin=Foam::min(Hmin,val);
		  //if (Pstream::master())
		  //  printf("proc=%d: iCell=%d/%d, H=%f\n",Pstream::myProcNo(),iC, mesh.nCells(),val);
		  if (val<0.0)
		    continue;
		  else if (val>dxC)
		    sum+=dxC*dx*dx;
		  else
		    sum+=val*dx*dx;
		}
	    }
      //Info<<"iCell="<<iC<<"/"<<mesh.nCells()<<", Hmax="<<Hmax+Pmin.z()<<",  Hmin="<<Hmin+Pmin.z()<<endl;
      if(sum/V>1.0)
	{
	  //Info<<"Something wrong: alpha="<<sum/V<<endl;
	  alpha1[iCell]=1.0;
	}
      else
	alpha1[iCell]=sum/V;
    }

  alpha1.correctBoundaryConditions();
}


///////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[])
{
//    #include "addRegionOption.H"
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createFields.H"
 
    #include "createTimeControls.H"
    #include "readTimeControls.H"
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"

  /*ijkZone ijkMesh(mesh);
  const label NCells=ijkMesh.Nx();
  const vector Pmin=ijkMesh.Pmin();  
  const vector Pmax=ijkMesh.Pmax();
  const scalar dx=ijkMesh.dx();*/

  scalar recTimeAlpha = 0;
  scalar recTimeNormal = 0;
  scalar recTimeCurvature = 0;
  //label nMax=2;
  //boolList interfaceCells(mesh.nCells(),false);
  //boolList atBoundary(mesh.nCells(),false);
  //ijkMesh.markBoundaryCells(atBoundary,2);

  //uniformStencil stencil(mesh,ijkMesh,nMax);
  
  word reconScheme=fvSolutionDict.get<word>("reconstructionScheme");
  autoPtr<reconstructionSchemes> surf =
    reconstructionSchemes::New(alpha1,phi,U,fvSolutionDict);

  word curvModel=transportProperties.get<word>("curvatureModel");
  interfaceForces surfForces(alpha1,U,transportProperties);



  while (runTime.run())
    {
      runTime++;
      	  
      for(int i=0;i<5;i++)
	{
	  Info<<"Setting alpha field..."<<endl;
	  Foam::autoPtr<waveFunction> waveField(Foam::waveFunction::New
				       (
					setAlphaFieldDict.get<word>("type"),
					setAlphaFieldDict
					));
      
	  const scalar interfaceTol=setAlphaFieldDict.getOrDefault<scalar>("interfaceTol",0.000001);      
	  setAlpha(mesh,waveField,alpha1);	  
	  recTimeAlpha= mesh.time().cpuTimeIncrement();
	  Info<<"recTimeAlpha="<<recTimeAlpha<<endl;		
	  
	  Info<<"Calculating normals..."<<endl;
	  surf->reconstruct();
	  const volVectorField& faceNormal = surf->normal();
	  recTimeNormal= mesh.time().cpuTimeIncrement();
	  Info<<"recTimeNormal="<<recTimeNormal<<endl;		

	  //calculate surface tension
	  Info<<"Calculating curvatures..."<<endl;
	  surfForces.correct();
	  surfForces.surfaceTensionForce();
	  recTimeCurvature= mesh.time().cpuTimeIncrement();
	  Info<<"recTimeCurvature="<<recTimeCurvature<<endl;		
	  volScalarField curvCalc = mesh.lookupObjectRef<volScalarField>("K_");
	}
      
      runTime.write();
    }

  Info<< "End\n" << endl;
  return 0;
}


// ************************************************************************* //
