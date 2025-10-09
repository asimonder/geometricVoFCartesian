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
    Foam::ijkZone

Description

    Developed by Asim Onder (2021-2025)

SourceFiles
    ijkZone.C

\*---------------------------------------------------------------------------*/

#ifndef ijkZone_H
#define ijkZone_H

#include "fvMesh.H"
#include "globalIndex.H"
#include "volFields.H"
#include "surfaceMesh.H"
#include "volMesh.H"

#include "IOobject.H"
#include "IOdictionary.H"
#include "MeshObject.H"
#include "IFstream.H"
#include "OFstream.H"
#include "cellSet.H"
#include "regIOobject.H" 


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class ijkZone Declaration
\*---------------------------------------------------------------------------*/

class ijkZone
:
    public regIOobject 
{

  private:

    // Private data
  //- Reference to mesh
  
    const fvMesh& mesh_;
    point Pmin_;
    point Pmax_;
    point PminL_;
    point PmaxL_;
  // Local window integer offsets
    label minIL_, minJL_, minKL_;
  label maxIL_, maxJL_, maxKL_;
  scalar dx_;scalar dy_;scalar dz_;
    label Nx_;label Ny_;label Nz_;
    label NxL_;label NyL_;label NzL_;
    scalar Lx_;scalar Ly_;scalar Lz_;
    labelList globalIds_;
    labelList neighbourProcs_;
    const globalIndex globalNumbering_;
    const volVectorField& C_;
    Vector<bool> isEmpty_;
    Vector<bool> isCyclic_;
    bool symXIn_;
    bool symXOut_;
    bool symYIn_;
    bool symYOut_;
    bool symZIn_;
    bool symZOut_;

  //template<class Type>
  //void collectDataFromOtherProcs(const List<labelHashSet>& commNeeded, Map<Type>& fieldData);

  //- Construct from fvMesh
    explicit ijkZone(const fvMesh&);
  
  //- No copy construct
    ijkZone(const ijkZone&) = delete;
  
    //- No copy assignment
    void operator=(const ijkZone&) = delete;

  public:
  // Constructors
  // Declare name of the class and its debug switch
    TypeName("ijkZone");
  
  // Selectors
  static ijkZone& New(const fvMesh&);
  
  //- Destructor
    virtual ~ijkZone() = default;      


  // Inherited functions
   bool writeData(Ostream&) const override;

  
    // Member functions
  
    void markBoundaryCells(boolList& atBoundary,label nOffsetCells);
    
    const labelList& globalIds() const
    {
      return globalIds_;
    } 

    void resetGlobalIds()
    {
      return globalIds_.clear();
    } 

    const Vector<bool> isEmpty() const
    {
      return isEmpty_;
    } 

    const Vector<bool> isCyclic() const
    {
      return isCyclic_;
    } 

    
    const globalIndex& globalNumbering() const
    {
      return globalNumbering_;
    } 

    template<typename Type>
    void getZoneField(const boolList& zone,
		      Map<Type>& phiMap,
		      const GeometricField<Type, fvPatchField, volMesh>& phi,
		      Vector<label> stencilSize);

    point Pmin() const {return Pmin_;}
    point Pmax() const {return Pmax_;}
    scalar dx() const {return dx_;}
    scalar dy() const {return dy_;}
    scalar dz() const {return dz_;}
    label Nx() const {return Nx_;};
    label Ny() const {return Ny_;};
    label Nz() const {return Nz_;};
    scalar Lx() const {return Lx_;};
    scalar Ly() const {return Ly_;};
    scalar Lz() const {return Lz_;};
    bool symXIn() const {return symXIn_;};
    bool symXOut() const {return symXOut_;};
    bool symYIn() const {return symYIn_;};
    bool symYOut() const {return symYOut_;};
    bool symZIn() const {return symZIn_;};
    bool symZOut() const {return symZOut_;};
    //convert i,j,k addressing to 1D 
    label ijk1(Vector<label> ijk3) const {return ijk3.x()+Nx_*ijk3.y()+Nx_*Ny_*ijk3.z();}
    label ijk1(label i, label j, label k) const {return i+Nx_*j+Nx_*Ny_*k;}
    label ijk1L(label i, label j, label k) const {return i+NxL_*j+NxL_*NyL_*k;}

    Vector<label> ijk3(const label celli) const;
    Vector<label> ijk1To3(label ijk) const; // const {return i+Nx_*j+Nx_*Ny_*k;}
    
  //    label ijkGlobalToLocal(const label celli);
    label ijkGlobalToLocal(const label ijk) const;
    label ijkGlobalToLocal(const label i, const label j, const label k) const;

    label lookupGlobalId(label i, label j, label k) const;
  
  };

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "ijkZoneI.H"
#endif


#endif

// ************************************************************************* //
