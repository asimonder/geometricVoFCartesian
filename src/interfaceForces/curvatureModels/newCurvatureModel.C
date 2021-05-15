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

\*---------------------------------------------------------------------------*/

#include "curvatureModel.H"
#include "messageStream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::curvatureModel>
Foam::curvatureModel::New
(
    const dictionary& dict,
    const volScalarField& alpha1,
    const volVectorField& U
)
{
    word modelType("divGradAlpha");

    if (!dict.readIfPresent("curvatureModel", modelType))
    {
        Warning
            << "Entry '"
            << "curvatureModel" << "' not found in dictionary "
            << dict.name() << nl
            << "using default" << nl;
    }

    Info<< "Selecting curvatureModel: " << modelType << endl;

    auto cstrIter = componentsConstructorTablePtr_->cfind(modelType);

    if (!cstrIter.found())
    {
        FatalIOErrorInLookup
        (
            dict,
            "curvatureModels",
            modelType,
            *componentsConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<curvatureModel>(cstrIter()( dict,alpha1, U));

}



// ************************************************************************* //
