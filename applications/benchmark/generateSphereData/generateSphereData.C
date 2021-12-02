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
    generateSphereData

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


template<class T>
void setAlpha
(
    T& cutCell,
    const fvMesh& mesh,
    const dictionary& initAlphaFieldDict,
    volScalarField& alpha1
)
{

    forAll(alpha1,cellI)
    {
        label cellStatus = cutCell.calcSubCell(cellI,0.0);

        if(cellStatus == -1)
        {
            alpha1[cellI] = 1;
        }
        else if(cellStatus == 1)
        {
            alpha1[cellI] = 0;
        }
        else if(cellStatus == 0)
        {
            if(cellStatus == 0 && mag(cutCell.faceArea()) != 0)
            {
                alpha1[cellI]= max(min(cutCell.VolumeOfFluid(),1),0);
            }
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

    #include "createTimeControls.H"
    #include "readTimeControls.H"
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"

  
    word reconScheme=fvSolutionDict.get<word>("reconstructionScheme");
    autoPtr<reconstructionSchemes> surf =
      reconstructionSchemes::New(alpha1,phi,U,fvSolutionDict);

    word curvModel=transportProperties.get<word>("curvatureModel");
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
    label nMax=1;
    boolList interfaceCells(mesh.nCells(),false);
    uniformStencil stencil(mesh,ijkMesh,nMax);
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

      runTime++;
      
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
      const volVectorField& faceNormal = surf->normal();
      volScalarField RDF = mesh.lookupObjectRef<volScalarField>("RDF");
      Map<scalar> rdfIJK;
      ijkMesh.getZoneField(interfaceCells,rdfIJK,RDF,stencilSize);

      //calculate surface tension
      surfForces.correct();
      surfForces.surfaceTensionForce();

      volScalarField curvCalc = mesh.lookupObjectRef<volScalarField>("K_");

      List<scalar> A;
      List<scalar> Df;
      scalar curv=1./R;
      scalar curvMax=0;

      //fprintf(writeTraining,"#G.x(),G.y(),curv*dx,A[0],A[1],A[2],A[3],A[4],A[5],A[6],A[7],A[8],A[9],A[10],A[11],A[12],A[13],A[14],A[15],A[16],A[17],A[18],A[19],A[20],A[21],A[22],A[23],A[24],A[25],A[26]\n");
      //fprintf(writeNormals,"#n.x(),n.y(),n.z()\n");
      //fprintf(writeCurvatures,"#curvature*dx\n");
      
      forAll(interfaceCells,iCell)		
	{
	  if(interfaceCells[iCell])
	    {
	      stencil.setStencil(alphaIJK,ijkMesh.ijk3(iCell));
	      A=stencil.getStencil();
	      stencil.setStencil(rdfIJK,ijkMesh.ijk3(iCell));
	      Df=stencil.getStencil();
	      Df=Df/dx;
	      vector n = faceNormal[iCell];
	      vector G=func.grad(mesh.C()[iCell]);

	      fprintf(writeTraining,"%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e\n",G.x(),G.y(),curv*dx,A[0],A[1],A[2],A[3],A[4],A[5],A[6],A[7],A[8],A[9],A[10],A[11],A[12],A[13],A[14],A[15],A[16],A[17],A[18],A[19],A[20],A[21],A[22],A[23],A[24],A[25],A[26],Df[0],Df[1],Df[2],Df[3],Df[4],Df[5],Df[6],Df[7],Df[8],Df[9],Df[10],Df[11],Df[12],Df[13],Df[14],Df[15],Df[16],Df[17],Df[18],Df[19],Df[20],Df[21],Df[22],Df[23],Df[24],Df[25],Df[26]);

	      fprintf(writeNormals,"%.12e,%.12e,%.12e\n",n.x(),n.y(),n.z());
	      fprintf(writeCurvatures,"%.12e\n",curvCalc[iCell]*dx);

	      curvMax=max(curvMax,curv);
	    }
	}	  
      Info<<"curvMax*dx="<<curvMax*dx<<endl;

      runTime.write();
    }


    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
