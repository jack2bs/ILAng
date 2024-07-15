
#include "ilang/ppa-estimations/ppa_hardware_block.h"
#include "ilang/ppa-estimations/ppa_profile_base.h"
#include <climits>
#include <ilang/ppa-estimations/ppa_model_registrar.h>
#include <ilang/ppa-estimations/ppa_profile_const_example.h>
#include <memory>

namespace ilang 
{

typedef PPAProfile_Const_Example constProf;
typedef std::shared_ptr<PPAProfile_Const_Example> constProf_ptr;

void registerAllModels(PPA_Registrar & registrar)
{

    // Uses the Aladdin power models for 10ns cycle time

    constProf_ptr ander_32b { 
        std::make_shared<constProf>(0.06, 1.807256e-03, 6.111633e-04, 5.036996e+01, 32, -1, None, 1)
    };

    constProf_ptr shifter_32b {
        std::make_shared<constProf>(0.70, 4.150581e-02, 1.695100e-03, 2.496461e+02, 32, 2, AcceleratorWide, 2)
    };

    constProf_ptr adder_32b {
        std::make_shared<constProf>(1.75, 8.534078e-03, 2.380803e-03, 1.794430e+02, 32, 4, AcceleratorWide, 2)
    };

    constProf_ptr multiplier_32b {
        std::make_shared<constProf>(4.9, 8.796144e-01, 4.817683e-02, 4.595000e+03, 32, 2, AcceleratorWide, 2)
    };

    constProf_ptr registerHold_32b {
        std::make_shared<constProf>(0, 0, 7.358945e-05, 5.981433e+00, 32, INT_MAX, None, 1)
    };

    constProf_ptr registerRead_32b {
        std::make_shared<constProf>(0, 5.535310e-04, 0, 0, 32, INT_MAX, None, 1)
    };

    constProf_ptr registerWrite_32b {
        std::make_shared<constProf>(0, 5.535310e-04, 0, 0, 32, INT_MAX, true, None, 1)
    };

    registrar.registerProfile(ander_32b, bBitwise);
    registrar.registerProfile(shifter_32b, bShifter);
    registrar.registerProfile(adder_32b, bAddition);
    // registrar.registerProfile(multiplier_32b, bDivision);
    registrar.registerProfile(multiplier_32b, bMultiplication);
    // registrar.registerProfile(multiplier_32b, bRemainder);
    registrar.registerProfile(registerHold_32b, bRegisterHold);
    registrar.registerProfile(registerRead_32b, bRegisterRead);
    registrar.registerProfile(registerWrite_32b, bRegisterWrite);

    registrar.finalizeRegistrar();
}


}