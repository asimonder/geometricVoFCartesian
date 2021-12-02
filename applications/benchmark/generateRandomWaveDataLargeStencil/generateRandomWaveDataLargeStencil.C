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
#include "writeFile.H"
#include "symmTensor2D.H"
//#include "isoCutCell.H"

#include "interfaceForces.H"
#include "ijkZone.H"
#include "uniformStencil.H"
#include "randomWaveFieldImplicitFunction.H"
#include "reconstructedDistanceFunction.H"
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
    const dictionary& setAlphaFieldDict,
    volScalarField& alpha1
)
{

  Foam::autoPtr<implicitFunction> func(Foam::implicitFunction::New
				       (
					setAlphaFieldDict.get<word>("type"),
					setAlphaFieldDict
					));	    

  const cellList& cells=mesh.cells();
  const faceList& faces=mesh.faces();
  const  surfaceVectorField&  Cf = mesh.Cf();
  //label Nx=40;
  label Nx=40;
  label Ny=Nx;
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
      Info<<"iCell="<<iC<<"/"<<mesh.nCells()<<", Hmax="<<Hmax+Pmin.z()<<",  Hmin="<<Hmin+Pmin.z()<<endl;
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

  ijkZone ijkMesh(mesh);
  const vector Pmin=ijkMesh.Pmin();  
  const vector Pmax=ijkMesh.Pmax();
  const scalar dx=ijkMesh.dx();
  const scalar pi=constant::mathematical::pi;
  word setAlphaMethod = fvSolutionDict.get<word>("setAlphaMethod");
  if (setAlphaMethod != "cutCellImpFunc" && setAlphaMethod != "cutCellIso")
    {
      FatalError  << "valid choice are only cutCellImpFunc or cutCellIso"
		  << abort(FatalError);
    }

  scalar amplitude=setAlphaFieldDict.get<scalar>("amplitude");
  scalar mwl=setAlphaFieldDict.get<scalar>("waterLevel");
  scalar interfaceTol=setAlphaFieldDict.get<scalar>("interfaceTol");
  label nRows=setAlphaFieldDict.get<label>("nRows");
  //setting the waves
  label Nk=setAlphaFieldDict.get<label>("Nk");
  label Nwaves=(2*Nk+1)*(Nk+1);
  List<scalar>kx(Nwaves,0.0);
  List<scalar>ky(Nwaves,0.0);
  List<scalar> comps(Nwaves,0.0);
  List<scalar> phases(Nwaves,0.0);

  label iC=0;
  for(int i=-Nk;i<Nk+1;i++)
    {
    for(int j=0;j<Nk+1;j++)
      {
	kx[iC]=2.0*pi*i; //L=1 hard coded!
	ky[iC]=2.0*pi*j;
	iC+=1;
      }
    }
  setAlphaFieldDict.add<List<scalar>>("kx",kx);
  setAlphaFieldDict.add<List<scalar>>("ky",ky);
  setAlphaFieldDict.add<List<scalar>>("comps",comps);
  setAlphaFieldDict.add<List<scalar>>("phases",phases);
  
  label seed=0;
  if (Pstream::master())
    seed=(unsigned) time(0);
  Foam::reduce(seed, maxOp<label>());
  Info<<"seed="<<(unsigned) time(0)<<endl;
  Random rndCentre(seed); 
  

  scalar recTime = 0;
  label nIter=setAlphaFieldDict.get<label>("nIter");;
  label nMax=(nRows-1)/2;
  boolList interfaceCells(mesh.nCells(),false);
  uniformStencil stencil(mesh,ijkMesh,nMax);
  
  autoPtr<reconstructionSchemes> surf =
    reconstructionSchemes::New(alpha1,phi,U,fvSolutionDict);

  interfaceForces surfForces(alpha1,U,transportProperties);
  
  while (runTime.run())
    {
      runTime++;
	  
      forAll(comps,i)
	{
	  comps[i]=rndCentre.globalPosition<scalar>(0,1); 
	  phases[i]=rndCentre.globalPosition<scalar>(0,2*constant::mathematical::pi); 
	}

      scalar compsL2=0;
      forAll(comps,i)
	compsL2+=Foam::sqr(comps[i]);

      Foam::reduce(compsL2, sumOp<scalar>());
      comps=comps/Foam::sqrt(compsL2)*amplitude;
	  
      setAlphaFieldDict.set<List<scalar>>("comps",comps);
      setAlphaFieldDict.set<List<scalar>>("phases",phases);
	  
      Info<<"Setting alpha field..."<<endl;
      setAlpha(mesh,setAlphaFieldDict,alpha1);	  
      Foam::implicitFunctions::randomWaveFieldImplicitFunction waveField(amplitude,
									 mwl,
									 phases,
									 comps,
									 kx,
									 ky);

      recTime += mesh.time().cpuTimeIncrement();
      Info<<"recTime="<<recTime<<endl;
	  
      Info<<"Calculating stencils..."<<endl;
      //select interfacecells
      interfaceCells=false;
      label cInt=0;
      forAll (alpha1,iCell)
	{
	  if (alpha1[iCell]>interfaceTol and alpha1[iCell]<1-interfaceTol)
	    {
	      interfaceCells[iCell]=true;
	      cInt+=1;
	    }
	}
      Info<<"number of interface cells="<<cInt<<endl;
      Map<scalar> alphaIJK;
      Vector<label> stencilSize(nMax,nMax,nMax);
      ijkMesh.getZoneField(interfaceCells,alphaIJK,alpha1,stencilSize);
      recTime += mesh.time().cpuTimeIncrement();
      Info<<"recTime="<<recTime<<endl;
      
      Info<<"Calculating curvatures..."<<endl;
      //interface reconstruction
      surf->reconstruct();
      const volVectorField& faceNormal = surf->normal();
      const boolList& interfaceCells = surf->interfaceCell();
      
      //calculate surface tension
      surfForces.correct();
      surfForces.surfaceTensionForce(); // update curvature	  
      recTime += mesh.time().cpuTimeIncrement();
      volScalarField curvHF = mesh.lookupObjectRef<volScalarField>("K_");
      
      List<scalar> A;
      scalar curv=0;
      scalar curvMax=0;

      FILE *writeObj;
      string fName="datasetsLargeStencil/curvDataset_amp"+std::to_string(amplitude)+"_Nk"+std::to_string(Nk)+"_"+std::to_string(Pstream::myProcNo())+".dat";
      if (Pstream::parRun())
	{
	  string dictPath=runTime.path()/".."/fName;
	  writeObj = fopen (dictPath.c_str(),"a");
	}
      else
	{
	  string dictPath=runTime.path()/fName;
	  writeObj = fopen (dictPath.c_str(),"a");
	}

      List <scalar>Hf;
      Hf.setSize(3*3*(2*nMax+1),0.0);  
      forAll(interfaceCells,iCell)		
	{
	  if(interfaceCells[iCell])
	    {
	      stencil.setStencil(alphaIJK,ijkMesh.ijk3(iCell));
	      A=stencil.getStencil();
	      
	      vector n = faceNormal[iCell];	      
	      label normalDir=-1;
	      if (mag(n.x())>0.0 and !ijkMesh.isEmpty().x())
		normalDir=0;
	      if (mag(n.y())>mag(n.x()) and !ijkMesh.isEmpty().y())
		normalDir=1;
	      if (mag(n.z())>mag(n.y()) and mag(n.z())>mag(n.x()) and !ijkMesh.isEmpty().z())
		normalDir=2;
	      
	      for (int iN=-nMax;iN<nMax+1;iN++)//normal 
		{
		  for (int iT=-1;iT<2;iT++)//tangent
		    {
		      for (int iB=-1;iB<2;iB++)//bi-normal
			{
			  if (normalDir==0)
			    {
			      Hf[(iT+1)+3*(iB+1)+9*(iN+nMax)]=A[stencil.a3(iN,iT,iB)];
			    }
			  else if (normalDir==1)
			    {
			      Hf[(iT+1)+3*(iB+1)+9*(iN+nMax)]=A[stencil.a3(iT,iN,iB)];
			    }
			  else
			    {
			      Hf[(iT+1)+3*(iB+1)+9*(iN+nMax)]=A[stencil.a3(iT,iB,iN)];
			    }			    			   
			}
		    }
		}
	      curv=waveField.curvature(mesh.C()[iCell]);
	      vector G=waveField.grad(mesh.C()[iCell]);
	      
	      fprintf(writeObj,"%.12e,%.12e,%.12e,%d,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e\n",G.x(),G.y(),curv*dx,normalDir,Hf[0],Hf[1],Hf[2],Hf[3],Hf[4],Hf[5],Hf[6],Hf[7],Hf[8],Hf[9],Hf[10],Hf[11],Hf[12],Hf[13],Hf[14],Hf[15],Hf[16],Hf[17],Hf[18],Hf[19],Hf[20],Hf[21],Hf[22],Hf[23],Hf[24],Hf[25],Hf[26],Hf[27],Hf[28],Hf[29],Hf[30],Hf[31],Hf[32],Hf[33],Hf[34],Hf[35],Hf[36],Hf[37],Hf[38],Hf[39],Hf[40],Hf[41],Hf[42],Hf[43],Hf[44],Hf[45],Hf[46],Hf[47],Hf[48],Hf[49],Hf[50],Hf[51],Hf[52],Hf[53],Hf[54],Hf[55],Hf[56],Hf[57],Hf[58],Hf[59],Hf[60],Hf[61],Hf[62]);

	      curvMax=max(curvMax,curv);
	    }
	}
      
      fclose(writeObj);

      Info<<"curvMax*dx="<<curvMax*dx<<endl;
      recTime += mesh.time().cpuTimeIncrement();
      Info<<"recTime="<<recTime<<endl;		
      
    }
  runTime.write();

  Info<< "End\n" << endl;
  return 0;
}


// ************************************************************************* //
