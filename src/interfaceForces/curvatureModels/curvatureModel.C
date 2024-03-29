/*---------------------------------------------------------------------------*\
            Copyright (c) 2021, Asim Onder 
-------------------------------------------------------------------------------
License
    This file is part of the VoFLibrary source code library, which is an
	unofficial extension to OpenFOAM.
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

#include "curvatureModel.H"
#include "surfaceInterpolate.H"
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(curvatureModel, 0);
    defineRunTimeSelectionTable(curvatureModel, components);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::curvatureModel::curvatureModel
(
    const word& type,
    const dictionary& dict,
    const volScalarField& alpha1,
    const volVectorField& U
)
:
    alpha1_(alpha1),
    U_(U),
    K_
    (
        IOobject
        (
            "K_",
            alpha1_.time().timeName(),
            alpha1_.mesh(),
	    IOobject::NO_READ,
	    //IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        alpha1_.mesh(),
        dimensionedScalar("K", dimless/dimLength, 0.0)//,
        //"zeroGradient"
	//"cyclic"
     ),
    Kf_(fvc::interpolate(K_))
{

}
    
// * * * * * * * * * * * * * * Public Access Member Functions  * * * * * * * * * * * * * * //


//    void Foam::curvatureModel::K();
//{
//    notImplemented("bool Foam::curvatureModel::correctContactAngle(scalar t)");;/
//}



// ************************************************************************* //
