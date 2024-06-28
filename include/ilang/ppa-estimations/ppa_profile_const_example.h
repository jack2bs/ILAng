#ifndef ILANG_ILA_PPA_ESTIMATIONS_PPA_PROFILE_CONST_EXAMPLE_H__
#define ILANG_ILA_PPA_ESTIMATIONS_PPA_PROFILE_CONST_EXAMPLE_H__

#include "ppa_profile_base.h"

namespace ilang {

class PPAProfile_Const_Example : public PPAProfile_Base
{
private:
    double m_blockTime;
    double m_blockDynamicPower;
    double m_blockLeakagePower;
    double m_blockArea;
    int m_maximumBitwidth;

public:
    PPAProfile_Const_Example
    (
        double blockTime,
        double blockDynamicPower,
        double blockLeakagePower,
        double blockArea,
        int maximumBitwidth
    ) : 
    m_blockTime(blockTime),
    m_blockDynamicPower(blockDynamicPower),
    m_blockLeakagePower(blockLeakagePower),
    m_blockArea(blockArea),
    m_maximumBitwidth(maximumBitwidth)
    {};

    double getBlockTime() override { return m_blockTime; }
    double getBlockDynamicPower() override { return m_blockDynamicPower; }
    double getBlockLeakagePower() override { return m_blockLeakagePower; }
    double getBlockArea() override { return m_blockArea; }
    int getMaximumBitwidth() override { return m_maximumBitwidth; }

};

}

#endif//ILANG_ILA_PPA_ESTIMATIONS_PPA_PROFILE_CONST_EXAMPLE_H__
