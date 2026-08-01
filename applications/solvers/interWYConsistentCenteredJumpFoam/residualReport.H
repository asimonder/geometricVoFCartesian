// residual face acceleration left by the projection, as a velocity
// increment dt*gf, over all internal faces and over faces touching a
// mixed cell.  Diagnostic only.
if (jumpReport && (runTime.timeIndex() % jumpReportInterval == 0))
{
    const labelUList& own = mesh.owner();
    const labelUList& nei = mesh.neighbour();
    const scalarField& a = alpha1.primitiveField();
    const scalar aTol = 1e-6;

    scalar sAll = 0, sMix = 0, gMaxMix = 0;
    label nAll = 0, nMix = 0;

    forAll(own, facei)
    {
        const scalar dv = runTime.deltaTValue()*gf[facei];
        sAll += dv*dv;
        ++nAll;

        const bool mixed =
            (a[own[facei]] > aTol && a[own[facei]] < 1 - aTol)
         || (a[nei[facei]] > aTol && a[nei[facei]] < 1 - aTol);

        if (mixed)
        {
            sMix += dv*dv;
            gMaxMix = max(gMaxMix, mag(dv));
            ++nMix;
        }
    }

    reduce(sAll, sumOp<scalar>());
    reduce(sMix, sumOp<scalar>());
    reduce(nAll, sumOp<label>());
    reduce(nMix, sumOp<label>());
    reduce(gMaxMix, maxOp<scalar>());

    Info<< "jump: residual dt*gf  rms all "
        << (nAll ? Foam::sqrt(sAll/nAll) : 0.0)
        << "  rms mixed " << (nMix ? Foam::sqrt(sMix/nMix) : 0.0)
        << "  max mixed " << gMaxMix
        << "  on " << nMix << " faces, maxU " << gMax(mag(U)()) << endl;
}
