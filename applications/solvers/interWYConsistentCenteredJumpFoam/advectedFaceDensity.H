// Face density from the density the sweeps advect.
//
//     alpha_f   = phi_alpha/phi, clipped to [0, 1]
//     rho_f     = rho1 alpha_f + rho2 (1 - alpha_f)
//
// Where phi passes through zero the ratio is meaningless, so the value is
// blended onto the linear face interpolation of alpha with the weight
// w = phi^2/(phi^2 + phiRef^2), phiRef a fixed fraction of max|phi|.
{
    const surfaceScalarField& aPhi = advector.alphaPhi();

    const dimensionedScalar phiSmall("phiSmall", phi.dimensions(), VSMALL);
    const dimensionedScalar phiRef
    (
        "phiRef",
        phi.dimensions(),
        advectedPhiTol*(gMax(mag(phi.primitiveField())) + VSMALL)
    );

    surfaceScalarField alphaf
    (
        min(max(aPhi/stabilise(phi, phiSmall), scalar(0)), scalar(1))
    );

    const surfaceScalarField w(sqr(phi)/(sqr(phi) + sqr(phiRef)));
    alphaf = w*alphaf + (scalar(1) - w)*fvc::interpolate(alpha1);

    rhof = rho1*alphaf + rho2*(scalar(1) - alphaf);

    if (jumpReport && (runTime.timeIndex() % jumpReportInterval == 0))
    {
        Info<< "jump: advected face density, rho_f in ["
            << gMin(rhof.primitiveField()) << ", "
            << gMax(rhof.primitiveField()) << "] kg/m3, phiRef "
            << phiRef.value() << endl;
    }
}
