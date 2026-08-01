/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2021-2025 Asim Onder
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
    Structured stencil on uniform Cartesian grids. Stores a flat array A_
    indexed by (i,j,k) offsets relative to a seed cell, with methods for
    interface normal estimation, flipping, and transposing.

    Subclassed by amrStencil for AMR grids where cells are discovered
    by topological traversal instead of global IJK indexing.

    Developed by Asim Onder (2021-2025)

SourceFiles
    uniformStencil.C

\*---------------------------------------------------------------------------*/

#ifndef uniformStencil_H
#define uniformStencil_H

#include "fvMesh.H"
#include "ijkZone.H"
#include "volFields.H"
#include "surfaceMesh.H"
#include "volMesh.H"

#include "IOobject.H"
#include "MeshObject.H"
#include "IFstream.H"
#include "cellSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class uniformStencil Declaration
\*---------------------------------------------------------------------------*/

  class uniformStencil
  {
  private:

    //- Reference to mesh
    const fvMesh& mesh_;

    //- Pointer to ijkZone (owned by registry, not by us).
    //  nullptr when constructed via SubclassTag (amrStencil).
    //  Used only by setStencil(phiIJK, ijk).
    const ijkZone* ijkMeshPtr_;


  protected:

    // Protected data -- accessible by amrStencil

    //- Total stencil size
    int N_;

    //- 2D flag
    bool is2D_;

    //- Flat stencil array
    List<scalar> A_;

    //- Stencil half-widths
    int iMax_;
    int jMax_;
    int kMax_;

    //- Stencil dimensions
    int NxS_;
    int NyS_;
    int NzS_;

    //- Empty direction flags
    Vector<bool> isEmpty_;

    //- Grid spacing (used by normal calculations)
    scalar dx_;

    //- Protected tag constructor for subclasses that don't use ijkZone.
    //  Initializes stencil storage and dimensions only; sets ijkMeshPtr_
    //  to nullptr.
    struct SubclassTag {};
    uniformStencil
    (
        SubclassTag,
        const fvMesh& mesh,
        int iMax,
        int jMax,
        int kMax,
        const Vector<bool>& isEmpty,
        scalar dx
    );


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
    virtual ~uniformStencil() = default;      

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
