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
#include "reconstructedDistanceFunction.H"
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
    volScalarField& alpha1,
    //List<vector>& refNormals,
    //scalarList& refCurvatures,
    bool is2D,
    scalar pertAmp
)
{

  const cellList& cells=mesh.cells();
  const faceList& faces=mesh.faces();
  const  surfaceVectorField&  Cf = mesh.Cf();
  const  volVectorField&  C = mesh.C();
  label Nx=40;
  label Ny=(is2D?1:Nx);
  label iC=0;
  Info<<"cells.size()="<<cells.size()<<endl;

  
  label seed=0;
  if (Pstream::master())
    seed=(unsigned) time(0);
  Foam::reduce(seed, maxOp<label>());
  Info<<"seed="<<(unsigned) time(0)<<endl;
  Random rndAmp(seed);

  forAll(alpha1,iCell)
    {
      iC+=1;
      vector Pmin=vector::uniform(VGREAT);
      vector Pmax=vector::uniform(-VGREAT);
      bool isWater=0;
      bool isAir=0;
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
      /*scalar distFac=1002.1;
      if (false)//only works in moderate to fine resolutions
	{
	  scalar H1=func->value(Pmin);
	  scalar H2=func->value(Pmax);
	  scalar Hmax=Foam::max(H1,H2);
	  scalar Hmin=Foam::min(H1,H2);
	  if (Pmax.z()+distFac*dxC<Hmin && Pmin.z()+distFac*dxC<Hmin)
	    {
	      alpha1[iCell]=1.0;
	      continue;
	    }
	  if (Pmax.z()-distFac*dxC>Hmax && Pmin.z()-distFac*dxC>Hmax)
	    {
	      alpha1[iCell]=0.0;
	      continue;
	    }
	    }*/

      //printf("proc=%d: iCell=%d/%d, dx=%f, dxC=%f\n",Pstream::myProcNo(),iC, mesh.nCells(),dx,dxC);
      
      scalar V=(is2D?dxC*dxC:dxC*dxC*dxC);
      scalar A=(is2D?dx:dx*dx);
      vector P0=Pmin+vector(dx/2.,dx/2.,0);
      scalar sum=0;
      vector meanNormal(0,0,0);
      vector Hw(0,0,0);
      scalar meanCurv=0;
      int nInterface=0;
      scalar Hmax=-VGREAT;
      scalar Hmin=VGREAT;
      for (int i=0;i<Nx;i++)
	    {
	      for (int j=0;j<Ny;j++)
		{
		  vector Pc=P0+vector(dx*i,dx*j,0);
		  scalar val=func->value(Pc);
		  //Hmax=Foam::max(Hmax,val);
		  //Hmin=Foam::min(Hmin,val);
		  //if (Pstream::master())
		  //  printf("proc=%d: iCell=%d/%d, H=%f\n",Pstream::myProcNo(),iC, mesh.nCells(),val);
		  if (val<0.0)
		    continue;
		  else if (val>dxC)
		    sum+=dxC*A;
		  else
		    {
		      sum+=val*A;
		      //vector G=func->grad(Pc);
		      //meanNormal+=G/Foam::mag(G);
		      //meanCurv+=func->curvature(Pc);
		      nInterface+=1;
		    }
		}
	    }
      if(sum/V>1.0)
	{
	  //Info<<"Something wrong: alpha="<<sum/V<<endl;
	  alpha1[iCell]=1.0;
	}
      else
	{
	  alpha1[iCell]=sum/V;
	  //Info<<"iCell="<<iC<<"/"<<mesh.nCells()<<endl;
	}

      /*  if (nInterface>0)
	{
	  refNormals[iCell]=meanNormal/float(nInterface);
	  refCurvatures[iCell]=meanCurv/float(nInterface);
	  //refNormals[iCell]=func->grad(C[iCell])/float(nInterface); //Foam::mag(func->grad(C[iCell]));
	  //refCurvatures[iCell]=func->curvature(C[iCell]);
	  //Info<<"refNormal="<<refNormals[iCell]<<", centerNormal="<<func->grad(C[iCell])/Foam::mag(func->grad(C[iCell]))<<endl;
	  //Info<<"refCurvature="<<refCurvatures[iCell]<<", centerCurvature="<<func->curvature(C[iCell])<<endl;
	  }*/
      
      //Info<<"iCell="<<iC<<"/"<<mesh.nCells()<<", Hmax="<<Hmax+Pmin.z()<<",  Hmin="<<Hmin+Pmin.z()<<endl;
      if (nInterface>0 && pertAmp>0.0)
	{
	  alpha1[iCell]+=alpha1[iCell]*rndAmp.globalPosition<scalar>(-pertAmp,pertAmp);
	  if(alpha1[iCell]>1.0)
	    {
	      alpha1[iCell]=1.0;
	    }
	  if(alpha1[iCell]<0.0)
	    {
	      alpha1[iCell]=0.0;
	    }
	}

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
  const label NCells=ijkMesh.Nx();
  const vector Pmin=ijkMesh.Pmin();  
  const vector Pmax=ijkMesh.Pmax();
  const scalar dx=ijkMesh.dx();

  bool is2D=ijkMesh.isEmpty().y()? true: false;
  
  scalar recTime = 0;
  label nMax=2;
  boolList interfaceCells(mesh.nCells(),false);
  //boolList atBoundary(mesh.nCells(),false);
  //ijkMesh.markBoundaryCells(atBoundary,2);

  uniformStencil stencil(mesh,ijkMesh,nMax);
  
  word reconScheme=fvSolutionDict.get<word>("reconstructionScheme");
  autoPtr<reconstructionSchemes> surf =
    reconstructionSchemes::New(alpha1,phi,U,fvSolutionDict);

  scalar pertAmp=fvSolutionDict.getOrDefault<scalar>("perturbationAmplitude",0.0);
  word curvModel=transportProperties.get<word>("curvatureModel");
  interfaceForces surfForces(alpha1,U,transportProperties);

  List<vector> refNormals(mesh.nCells());
  scalarList refCurvatures(mesh.nCells(),0);

  //reconstructedDistanceFunction RDF(mesh);
  
  while (runTime.run())
    {
      runTime++;
      	  
      Info<<"Setting alpha field..."<<endl;
      //setAlpha(cutCell,mesh,setAlphaFieldDict,alpha1);
      Foam::autoPtr<waveFunction> waveField(Foam::waveFunction::New
				       (
					setAlphaFieldDict.get<word>("type"),
					setAlphaFieldDict
					));
      bool randomHeight=setAlphaFieldDict.getOrDefault<bool>("randomHeight",false);
      if (randomHeight)
	{
	  label seed=0;
	  if (Pstream::master())
	    seed=(unsigned) time(0);
	  Foam::reduce(seed, maxOp<label>());
	  Random rndWaterLevel(seed); 
	  scalar Hw=rndWaterLevel.globalPosition<scalar>(0,dx);
	  waveField->setWaterLevel(Hw);
	  Info<<"Hw="<<Hw<<", Hw/dx="<<Hw/dx<<endl;
	}
      const scalar interfaceTol=setAlphaFieldDict.getOrDefault<scalar>("interfaceTol",0.000001);      
      setAlpha(mesh,waveField,alpha1,is2D,pertAmp);
      recTime += mesh.time().cpuTimeIncrement();

      FILE *writeTraining;
      FILE *writeNormals;
      FILE *writeCurvatures;
      string strNc=std::to_string(NCells);
      //string strXFetch=std::to_string(XFetch);
      //string strU=std::to_string(UWind);
      string strName=waveField->name();
      string dictPath=runTime.path()/"trainingData_"+strName+"_N"+strNc+".dat";
      word mlModel=fvSolutionDict.getOrDefault<word>("mlModel","alpha");
      if (mlModel=="rdf")
	dictPath=runTime.path()/"trainingDataRDF_"+strName+"_N"+strNc+".dat";
      if (randomHeight)
	dictPath=runTime.path()/"trainingDataRandomHeight_"+strName+"_N"+strNc+".dat";

      writeTraining = fopen (dictPath.c_str(),"a");
      dictPath=runTime.path()/"normalData_"+strName+"_N"+strNc+"_"+reconScheme+".dat";
      writeNormals = fopen (dictPath.c_str(),"a");
      dictPath=runTime.path()/"curvatureData_"+strName+"_N"+strNc+"_"+curvModel+"_"+reconScheme+".dat";
      writeCurvatures = fopen (dictPath.c_str(),"a");

      Info<<"recTime="<<recTime<<endl;
      
      Info<<"Calculating stencils..."<<endl;
	  //select interfacecells
      interfaceCells=false;
      label cInt=0;
      forAll (alpha1,iCell)
	{
	  if (alpha1[iCell]>interfaceTol and alpha1[iCell]<1-interfaceTol) // and !atBoundary[iCell])
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
      volScalarField RDF = mesh.lookupObjectRef<volScalarField>("RDF");
      Map<scalar> rdfIJK;
      ijkMesh.getZoneField(interfaceCells,rdfIJK,RDF,stencilSize);

      //calculate surface tension
      surfForces.correct();
      surfForces.surfaceTensionForce();

      recTime += mesh.time().cpuTimeIncrement();
      volScalarField curvCalc = mesh.lookupObjectRef<volScalarField>("K_");     
      
      List<scalar> A;
      scalar curv=0;
      scalar curvMax=0;

      const volVectorField& C = mesh.C();
      
      forAll(interfaceCells,iCell)		
	{
	  if(interfaceCells[iCell])
	    {
	      if (mlModel=="rdf")
		{
		  stencil.setStencil(rdfIJK,ijkMesh.ijk3(iCell));
		  A=stencil.getStencil();
		  A=A/dx;
		}
	      else
		{
		  stencil.setStencil(alphaIJK,ijkMesh.ijk3(iCell));
		  A=stencil.getStencil();
		}
	      vector n = faceNormal[iCell];	      

	      vector X=C[iCell];
	      curv=waveField->curvature(X);
	      vector G=waveField->grad(X);
	      //curv=refCurvatures[iCell];
	      //vector G=refNormals[iCell]/Foam::mag(refNormals[iCell].z());
	      //Info<<"G.x()="<<G.x()<<", Gc.x()="<<Gc.x()<<endl;
	      symmTensor2D H=waveField->hessian(X);
	      //if (is2D)
	      //fprintf(writeTraining,"%.12e,%.12e,%.12e,%.12e,",G.x(),G.y(),curv*dx,dx);
	      //else
	      fprintf(writeTraining,"%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,",G.x(),G.y(),H.xx()*dx,H.yy()*dx,H.xy()*dx,curv*dx,dx);
	      for(int iA =0;iA<A.size()-1;iA++)
		  fprintf(writeTraining,"%.12e,",A[iA]);
	      fprintf(writeTraining,"%.12e\n",A[A.size()-1]);
	      n/=Foam::mag(n);
	      //G=refNormals[iCell]/Foam::mag(refNormals[iCell]);
	      G/=Foam::mag(G);
	      fprintf(writeNormals,"%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e\n",X.x(),X.y(),X.z(),n.x(),n.y(),n.z(),G.x(),G.y(),G.z());
	      fprintf(writeCurvatures,"%.12e,%.12e,%.12e,%.12e,%.12e\n",X.x(),X.y(),X.z(),curvCalc[iCell]*dx,curv*dx);

	      curvMax=max(curvMax,curv);
	    }
	}	  
      Info<<"curvMax*dx="<<curvMax*dx<<endl;
      recTime += mesh.time().cpuTimeIncrement();
      Info<<"recTime="<<recTime<<endl;		

      fclose(writeTraining);
      fclose(writeNormals);
      fclose(writeCurvatures);
	
      runTime.write();
    }

  Info<< "End\n" << endl;
  return 0;
}


// ************************************************************************* //
