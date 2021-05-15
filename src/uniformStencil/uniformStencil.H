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
    List<scalar> A_;
    const label nMax_;
    label N_;
    label NS_;
    bool is2D_;
      
  //- No copy construct
    uniformStencil(const uniformStencil&) = delete;
  
    //- No copy assignment
    void operator=(const uniformStencil&) = delete;

    label a2(const label i, const label j)
    {
      return (i+1)+NS_*(j+1);
    }
  
    label a3(const label i, const label j, const label k)
    {
      return (i+1)+NS_*(j+1)+NS_*NS_*(k+1);
    }

  public:
  // Constructors
  // Declare name of the class and its debug switch
    TypeName("uniformStencil");

  //- Construct from fvMesh
    explicit uniformStencil(const fvMesh&,const ijkZone& ijkMesh,const label nMax);

  // Selectors
  //static uniformStencil& New(const fvMesh&);

  //- Destructor
    virtual ~uniformStencil() = default;      

    void setStencil(const Map<scalar>& phiIJK,const Vector<label>& ijk);

    const List<scalar> getStencil(){return A_;};

    vector calcYoungNormal();
    vector calcCCDNormal();
    
  };

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//#ifdef NoRepository
//#include "uniformStencilI.H"
//#endif


#endif

// ************************************************************************* //
