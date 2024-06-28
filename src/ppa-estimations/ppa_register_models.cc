
#include "ilang/ppa-estimations/ppa_hardware_block.h"
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
    constProf_ptr nohardware {
        std::make_shared<constProf>(0, 0, 0, 0, INT_MAX)
    };

    constProf_ptr ander_32b { 
        std::make_shared<constProf>(3, 1611.753, 493.462, 9.437, 32) 
    };

    constProf_ptr shifter_32b {
        std::make_shared<constProf>(74, 17331.483, 2270.682, 64.684, 32)
    };

    constProf_ptr adder_32b {
        std::make_shared<constProf>(302, 6711.251, 942.731, 36.913, 32)
    };

    constProf_ptr divider_32b {
        std::make_shared<constProf>(5356, 189528.282, 26600.057, 804.422, 32)
    };

    constProf_ptr multiplier_32b {
         std::make_shared<constProf>(355, 228240.874, 17879.126, 509.903, 32)
    };
    constProf_ptr remainder_32b { 
        std::make_shared<constProf>(5612, 203004.908, 27450.345, 833.176, 32)
    };

    registrar.registerProfile(nohardware, bNoHardware);
    registrar.registerProfile(ander_32b, bBitwise);
    registrar.registerProfile(shifter_32b, bShifter);
    registrar.registerProfile(adder_32b, bAddition);
    registrar.registerProfile(divider_32b, bDivision);
    registrar.registerProfile(multiplier_32b, bMultiplication);
    registrar.registerProfile(remainder_32b, bRemainder);


    registrar.finalizeRegistrar();
}


}