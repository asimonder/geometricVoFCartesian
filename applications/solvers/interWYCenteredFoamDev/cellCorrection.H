/*---------------------------------------------------------------------------*\
    Part of interWYCenteredFoam (see the main file for the full
    description and license). Extracted for the conventional
    fragmented solver layout.
\*---------------------------------------------------------------------------*/

        // --- 5. averaged cell correction (approximate projection) -----
        const surfaceScalarField gf
        (
            af - (alphaf/dtDim)*fvc::snGrad(p_rgh)
        );
        faceToCellAvg(gf, gC);
        if (zonal)
        {
            // outside the uniform zone the axis-average is not
            // defined; use the standard consistent reconstruction
            // there (safe: density is constant outside, so the face
            // content carries no interfacial balance)
            const volVectorField gRec
            (
                fvc::reconstruct(gf*mesh.magSf())
            );
            const boolList& zc = splitPtr->inZone();
            vectorField& gcf = gC.primitiveFieldRef();
            forAll(gcf, celli)
            {
                if (!zc[celli]) gcf[celli] = gRec[celli];
            }
            gC.correctBoundaryConditions();
        }

        U += dtDim*gC;
        U.correctBoundaryConditions();

