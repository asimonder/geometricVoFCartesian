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

Class
    Foam::fv::momentumConsistentCrankNicolson

Description
    Hybrid temporal scheme that uses Euler at the interface (detected via 
    snGrad(alpha)) and Crank-Nicolson in bulk phases for consistency with
    isoAdvector in large density ratio flows.
    
    Usage:
    \verbatim
    ddtSchemes
    {
        default  momentumConsistentCrankNicolson 0.9 0.01;
        // 0.9 is bulk Crank-Nicolson coefficient
        // 0.01 is interface threshold for |grad(alpha)|
    }
    \endverbatim

\*---------------------------------------------------------------------------*/

#ifndef momentumConsistentCrankNicolson_H
#define momentumConsistentCrankNicolson_H

#include "ddtScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                  Class momentumConsistentCrankNicolson Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class momentumConsistentCrankNicolson
:
    public fv::ddtScheme<Type>
{
    // Private Data

        //- Store ddt0 fields on the registry
        template<class GeoField>
        class DDt0Field
        :
            public GeoField
        {
            label startTimeIndex_;

        public:

            //- Construct from file
            DDt0Field
            (
                const IOobject& io,
                const fvMesh& mesh
            );

            //- Construct value-initialised
            DDt0Field
            (
                const IOobject& io,
                const fvMesh& mesh,
                const typename GeoField::value_type& value,
                const dimensionSet& dims
            );

            //- Return start time index
            label startTimeIndex() const noexcept { return startTimeIndex_; }

            //- Cast to underlying GeoField
            GeoField& operator()();

            //- Assignment
            void operator=(const GeoField& gf);
        };

        //- Reference to alpha field
        const volScalarField& alpha_;
        
        //- Bulk Crank-Nicolson coefficient
        scalar bulkCnCoeff_;
        
        //- Interface threshold for |grad(alpha)|
        scalar interfaceThreshold_;
        
        //- Spatially varying offset coefficient field
        mutable scalarField ocCoeff_;


    // Private Member Functions

        //- No copy construct
        momentumConsistentCrankNicolson(const momentumConsistentCrankNicolson&) = delete;

        //- No copy assignment
        void operator=(const momentumConsistentCrankNicolson&) = delete;
        
        //- Update spatially varying coefficient
        void updateOcCoeff() const;
        
        //- Get or create ddt0 field
        template<class GeoField>
        DDt0Field<GeoField>& ddt0_
        (
            const word& name,
            const dimensionSet& dims
        );

        //- Check if ddt0 needs evaluation
        template<class GeoField>
        bool evaluate(DDt0Field<GeoField>& ddt0);

        //- Return coefficient (spatially varying)
        template<class GeoField>
        scalarField coef_(const DDt0Field<GeoField>&) const;

        //- Return old time coefficient (spatially varying)
        template<class GeoField>
        scalarField coef0_(const DDt0Field<GeoField>&) const;

        //- Return reciprocal time coefficient
        template<class GeoField>
        tmp<scalarField> rDtCoef_(const DDt0Field<GeoField>&) const;

        //- Return reciprocal old time coefficient
        template<class GeoField>
        tmp<scalarField> rDtCoef0_(const DDt0Field<GeoField>&) const;

        //- Return ddt0 multiplied by off-centering coefficient
        template<class GeoField>
        tmp<GeoField> offCentre_(const GeoField& ddt0) const;


public:

    //- Runtime type information
    TypeName("momentumConsistentCrankNicolson");


    // Constructors

        //- Construct from mesh and Istream
        momentumConsistentCrankNicolson(const fvMesh& mesh, Istream& is);


    // Member Functions

        //- Return mesh reference
        const fvMesh& mesh() const
        {
            return fv::ddtScheme<Type>::mesh();
        }

        tmp<GeometricField<Type, fvPatchField, volMesh>> fvcDdt
        (
            const dimensioned<Type>&
        );

        tmp<GeometricField<Type, fvPatchField, volMesh>> fvcDdt
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        tmp<GeometricField<Type, fvPatchField, volMesh>> fvcDdt
        (
            const dimensionedScalar&,
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        tmp<GeometricField<Type, fvPatchField, volMesh>> fvcDdt
        (
            const volScalarField&,
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        tmp<GeometricField<Type, fvPatchField, volMesh>> fvcDdt
        (
            const volScalarField& alpha,
            const volScalarField& rho,
            const GeometricField<Type, fvPatchField, volMesh>& psi
        );

        tmp<fvMatrix<Type>> fvmDdt
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        tmp<fvMatrix<Type>> fvmDdt
        (
            const dimensionedScalar&,
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        tmp<fvMatrix<Type>> fvmDdt
        (
            const volScalarField&,
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        tmp<fvMatrix<Type>> fvmDdt
        (
            const volScalarField& alpha,
            const volScalarField& rho,
            const GeometricField<Type, fvPatchField, volMesh>& psi
        );

        typedef typename ddtScheme<Type>::fluxFieldType fluxFieldType;

        tmp<fluxFieldType> fvcDdtUfCorr
        (
            const GeometricField<Type, fvPatchField, volMesh>& U,
            const GeometricField<Type, fvsPatchField, surfaceMesh>& Uf
        );

        tmp<fluxFieldType> fvcDdtPhiCorr
        (
            const GeometricField<Type, fvPatchField, volMesh>& U,
            const fluxFieldType& phi
        );

        tmp<fluxFieldType> fvcDdtUfCorr
        (
            const volScalarField& rho,
            const GeometricField<Type, fvPatchField, volMesh>& U,
            const GeometricField<Type, fvsPatchField, surfaceMesh>& Uf
        );

        tmp<fluxFieldType> fvcDdtPhiCorr
        (
            const volScalarField& rho,
            const GeometricField<Type, fvPatchField, volMesh>& U,
            const fluxFieldType& phi
        );

        tmp<surfaceScalarField> meshPhi
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "momentumConsistentCrankNicolson.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
