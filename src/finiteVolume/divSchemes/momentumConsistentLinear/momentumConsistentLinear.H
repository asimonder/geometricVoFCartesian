/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2025 Asim Onder
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

#include "upwind.H"
#include "volFields.H"
#include "immiscibleIncompressibleTwoPhaseGenericMixture.H"
//#include "incompressibleInterPhaseTransportModel.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifndef momentumConsistentLinear_H
#define momentumConsistentLinear_H

namespace Foam
{

template<class Type>
class momentumConsistentLinear
:
    public upwind<Type>
{
    // Private Data
        const volScalarField& alpha_;
        scalar interfaceThreshold_;

    // Private Member Functions
        momentumConsistentLinear(const momentumConsistentLinear&) = delete;
        void operator=(const momentumConsistentLinear&) = delete;

public:

    //- Runtime type information
    TypeName("momentumConsistentLinear");

    // Constructors
        momentumConsistentLinear
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            Istream& schemeData
        )
	  :
            upwind<Type>(mesh, faceFlux, schemeData),
	    alpha_
	    (
	     [&mesh]() -> const volScalarField&
	     {
	       // Try common alpha field names
	       if (mesh.foundObject<volScalarField>("alpha.phase1"))
		 return mesh.lookupObject<volScalarField>("alpha.phase1");
	       else if (mesh.foundObject<volScalarField>("alpha.water"))
		 return mesh.lookupObject<volScalarField>("alpha.water");
	       else if (mesh.foundObject<volScalarField>("alpha.liquid"))
		 return mesh.lookupObject<volScalarField>("alpha.liquid");
	       else if (mesh.foundObject<volScalarField>("alpha1"))
		 return mesh.lookupObject<volScalarField>("alpha1");
	       else
		 {
		   FatalErrorInFunction
		     << "Cannot find alpha field (tried alpha.phase1, alpha.water, alpha1)"
		     << abort(FatalError);
		   return mesh.lookupObject<volScalarField>("alpha.phase1"); // Dummy
		 }
	     }()
	     ),
            interfaceThreshold_(readScalar(schemeData))
        {}

     momentumConsistentLinear
      (
        const fvMesh& mesh,
        Istream& schemeData
       )
	:
        upwind<Type>(mesh, schemeData),
	alpha_
	(
	 [&mesh]() -> const volScalarField&
	     {
	       // Try common alpha field names
	       if (mesh.foundObject<volScalarField>("alpha.phase1"))
		 return mesh.lookupObject<volScalarField>("alpha.phase1");
	       else if (mesh.foundObject<volScalarField>("alpha.water"))
		 return mesh.lookupObject<volScalarField>("alpha.water");
	       else if (mesh.foundObject<volScalarField>("alpha.liquid"))
		 return mesh.lookupObject<volScalarField>("alpha.liquid");
	       else if (mesh.foundObject<volScalarField>("alpha1"))
		 return mesh.lookupObject<volScalarField>("alpha1");
	       else
		 {
		   FatalErrorInFunction
		     << "Cannot find alpha field (tried alpha.phase1, alpha.water, alpha1)"
		     << abort(FatalError);
		   return mesh.lookupObject<volScalarField>("alpha.phase1"); // Dummy
		 }
	     }()
	 ),
        interfaceThreshold_(readScalar(schemeData))
    {}
  
    // Member Functions
        virtual tmp<surfaceScalarField> weights
        (
            const GeometricField<Type, fvPatchField, volMesh>& vf
        ) const;
};

} // End namespace Foam

#endif
