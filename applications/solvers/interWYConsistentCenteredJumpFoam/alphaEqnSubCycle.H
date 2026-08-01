if (pimple.nCorrPIMPLE() > 1 || nAlphaSubCycles > 1)
{
    FatalErrorInFunction
        << "interWYConsistentCenteredFoam admits no outer correctors and"
        << " no alpha sub-cycles: the consistency of the cycle is"
        << " structural. Set nOuterCorrectors 1 and nAlphaSubCycles 1."
        << exit(FatalError);
}

#include "alphaEqn.H"
