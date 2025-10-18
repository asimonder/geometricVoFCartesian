#include "incompressibleInterPhaseTransportModel.H"
#include "immiscibleIncompressibleTwoPhaseGenericMixture.H"

namespace Foam
{
    typedef incompressibleInterPhaseTransportModel
    <
        immiscibleIncompressibleTwoPhaseGenericMixture
    > immiscibleIncompressibleTwoPhaseGenericMixtureIncompressibleInterPhaseTransportModel;

    defineTemplateTypeNameAndDebugWithName
    (
        immiscibleIncompressibleTwoPhaseGenericMixtureIncompressibleInterPhaseTransportModel,
        "immiscibleIncompressibleTwoPhaseGenericMixture<incompressibleInterPhaseTransportModel>",
        0
    );
}
