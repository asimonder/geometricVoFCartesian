// Count upwind faces
label nUpwindFaces = 0;
label nTotalFaces = mesh.nInternalFaces();

forAll(interfaceMask.primitiveField(), facei)
{
    if (interfaceMask.primitiveField()[facei] > 0.5)
    {
        nUpwindFaces++;
    }
}

reduce(nUpwindFaces, sumOp<label>());
reduce(nTotalFaces, sumOp<label>());

scalar upwindRatio = scalar(nUpwindFaces) / scalar(nTotalFaces);

Info<< "Upwind faces: " << nUpwindFaces << " / " << nTotalFaces 
    << " (" << 100*upwindRatio << "%)" << endl;
