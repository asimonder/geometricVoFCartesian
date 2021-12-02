IOdictionary setAlphaFieldDict
(
 IOobject
 (
  "setAlphaFieldDict",
  mesh.time().system(),
  mesh,
  IOobject::MUST_READ,
  IOobject::AUTO_WRITE
  )
 );


Info<< "Reading field U\n" << endl;
volVectorField U
(
 IOobject
 (
  "U",
  runTime.timeName(),
  mesh,
  IOobject::MUST_READ,
  IOobject::NO_WRITE
  ),
 mesh
 );

#include "createPhi.H"

Info<< "Reading field alpha1\n" << endl;
volScalarField alpha1
(
 IOobject
 (
  "alpha1",
  runTime.timeName(),
  mesh,
  IOobject::MUST_READ,
  IOobject::AUTO_WRITE
  ),
 mesh
 );

IOdictionary transportProperties
(
 IOobject
 (
  "transportProperties",
  U.time().constant(),
  U.db(),
  IOobject::MUST_READ_IF_MODIFIED,
  IOobject::NO_WRITE
  )
 );


IOdictionary fvSolutionDict
(
 IOobject
 (
  "fvSolution",
  alpha1.time().system(),
  alpha1.db(),
  IOobject::MUST_READ_IF_MODIFIED,
  IOobject::NO_WRITE
  )
 );
