/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2021 Asim Onder
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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

Class
    Foam::uniformStencil

Description
    Developed by Asim Onder (2021)

SourceFiles
    uniformStencil.C

\*---------------------------------------------------------------------------*/

#ifndef uniformStencil_H
#define uniformStencil_H

#include "fvMesh.H"
#include "ijkZone.H"
//#include "globalIndex.H"
#include "volFields.H"
#include "surfaceMesh.H"
#include "volMesh.H"

#include "IOobject.H"
#include "MeshObject.H"
#include "IFstream.H"
#include "cellSet.H"
//#include "DynamicList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class uniformStencil Declaration
\*---------------------------------------------------------------------------*/

  class uniformStencil
  {
  private:

    // Private data
  //- Reference to mesh
  
    const fvMesh& mesh_;
    const ijkZone& ijkMesh_;
    int N_;
    bool is2D_;
    List<scalar> A_;
    int iMax_;
    int jMax_;
    int kMax_;
    int NxS_;
    int NyS_;
    int NzS_;

  //- No copy construct
    //uniformStencil(const uniformStencil&) = delete;
  
    //- No copy assignment
    //void operator=(const uniformStencil&) = delete;

  public:
  // Constructors
    //uniformStencil();
    uniformStencil() = default;
    
  //- Construct from fvMesh
    uniformStencil(const fvMesh&,const ijkZone& ijkMesh,int nMax);
    uniformStencil(const fvMesh&,const ijkZone& ijkMesh,int iMax,int jMax,int kMax);

  // Selectors
  //static uniformStencil& New(const fvMesh&);

  //- Destructor
    // virtual ~uniformStencil() = default;      

    void setStencil(const Map<scalar>& phiIJK,const Vector<label>& ijk);

    const List<scalar> getStencil(){return A_;};

    //scalar estimateSignK();
    //scalar estimateSignK(vector n);
    vector calcYoungNormal();
    vector calcCCDNormal();
    //scalar calcCurvature(const Map<scalar>& phiIJK,const Vector<label>& ijk,label normalDir);
    //scalar calcCurvature(label normalDir);

    void flip(const int axis);
    void transpose(const bool xaxis,const bool yaxis,const bool zaxis);
    List<scalar> getSmallStencil(const int iMax,const int jMax,const int kMax,std::string format);

    label a2(const label i, const label j)
    {
      //hardcoded to z!
      return (i+iMax_)+NxS_*(j+kMax_);
    }
  
    label a3(const label i, const label j, const label k)
    {
      //return (k+nMax_)+NxS_*(j+nMax_)+NxS_*NyS_*(i+nMax_);
      return (i+iMax_)+NxS_*(j+jMax_)+NxS_*NyS_*(k+kMax_);
    }

    //   std::vector<double> getSmallStencil(int iMin, int jMin, int kMin);
  };

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//#ifdef NoRepository
//#include "uniformStencilI.H"
//#endif


#endif

// ************************************************************************* //
