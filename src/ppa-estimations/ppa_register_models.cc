
#include "ilang/ppa-estimations/ppa_hardware_block.h"
#include <ilang/ppa-estimations/ppa_model_registrar.h>
#include <ilang/ppa-estimations/ppa_profile_const_example.h>
#include <memory>
#include <climits>
#include <ilang/ppa-estimations/ppa_callbacks.h>

namespace ilang 
{

typedef PPAProfile_Const_Example constProf;
typedef std::shared_ptr<PPAProfile_Const_Example> constProf_ptr;

void registerAllModels(PPA_Registrar & registrar)
{

    constProf_ptr ander_32b { 
        std::make_shared<constProf>(3, 1611.753, 493.462, 9.437, 32, -1, false) 
    };

    constProf_ptr shifter_32b {
        std::make_shared<constProf>(74, 17331.483, 2270.682, 64.684, 32, 400, true)
    };

    constProf_ptr adder_32b {
        std::make_shared<constProf>(302, 6711.251, 942.731, 36.913, 32, 400, true)
    };

    constProf_ptr divider_32b {
        std::make_shared<constProf>(5356, 189528.282, 26600.057, 804.422, 32, 1, true)
    };

    constProf_ptr multiplier_32b {
         std::make_shared<constProf>(355, 228240.874, 17879.126, 509.903, 32, 23, true)
    };

    constProf_ptr remainder_32b { 
        std::make_shared<constProf>(5612, 203004.908, 27450.345, 833.176, 32, 1, true)
    };

    class ProfRegisterHold : public PPAProfile_Base
    {
    private:
        const int m_bitwidth;

    public:
        ProfRegisterHold(int bitwidth) : m_bitwidth(bitwidth){}

        double getBlockTime() override { return 0.0; }

        virtual double getBlockDynamicPower() override { return 0.0; }

        virtual double getBlockLeakagePower() override 
            { return 47.339 * m_bitwidth; }

        virtual double getBlockArea() override 
            { return 1.278 * m_bitwidth; }

        virtual int getMaximumBitwidth() override { return m_bitwidth; }

        virtual int getMaximumInstances() override { return INT_MAX; }

        virtual bool getIsReusable() override { return false; }
    };



    class ProfRegisterRW : public PPAProfile_Base
    {
    private:
        const double m_dynPowerConst;
        const int m_bitwidth;

    public:
        ProfRegisterRW(double dynPowerConst, int bitwidth) 
            : m_dynPowerConst(dynPowerConst), m_bitwidth(bitwidth){}

        double getBlockTime() override { return 0.0; }

        virtual double getBlockDynamicPower() override 
        { 
            double switches = PPA_Callbacks::getSwitchingActivity();
            return switches * m_dynPowerConst;
        }

        virtual double getBlockLeakagePower() override { return 0.0; }

        virtual double getBlockArea() override { return 0.0; }

        virtual int getMaximumBitwidth() override { return m_bitwidth; }

        virtual int getMaximumInstances() override { return INT_MAX; }

        virtual bool getIsReusable() override { return false; }
    };

    for (int i = 1; i <= 256; i++)
    {
        std::shared_ptr<ProfRegisterRW> register_read =
            std::make_shared<ProfRegisterRW>((404.466 + 19.181) * i, i);
        std::shared_ptr<ProfRegisterRW> register_write =
            std::make_shared<ProfRegisterRW>((404.466 + 19.181) * i, i);
        std::shared_ptr<ProfRegisterHold> register_hold =
            std::make_shared<ProfRegisterHold>(i);

        registrar.registerProfile(register_read, bRegisterRead);
        registrar.registerProfile(register_write, bRegisterWrite);
        registrar.registerProfile(register_hold, bRegisterHold);
    }



    registrar.registerProfile(ander_32b, bBitwise);
    registrar.registerProfile(shifter_32b, bShifter);
    registrar.registerProfile(adder_32b, bAddition);
    registrar.registerProfile(divider_32b, bDivision);
    registrar.registerProfile(multiplier_32b, bMultiplication);
    registrar.registerProfile(remainder_32b, bRemainder);


    registrar.finalizeRegistrar();
}


}