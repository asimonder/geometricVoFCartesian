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

\*---------------------------------------------------------------------------*/

#include "ijkZone.H"


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
template<typename Type>
void Foam::ijkZone::getZoneField(const Foam::boolList& zone,
				Foam::Map<Type>& phiIJK,
				 const Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh>& phi,
				 Foam::Vector<label> stencilSize)
{

  List<labelHashSet> commNeeded(Pstream::nProcs());
  List<labelHashSet> send(Pstream::nProcs());
  //const volVectorField& C = mesh_.C();
  label iMax=stencilSize.x();
  label jMax=stencilSize.y();
  label kMax=stencilSize.z();
  //bool is2D=false;

  if (isEmpty_.x())
    iMax=0;
  if (isEmpty_.y())
    jMax=0;
  if (isEmpty_.z())
    kMax=0;
  
#include "ijkZone/zoneTransfers.H"
  
}


// ************************************************************************* //
