/// \file
/// 


#ifndef ILANG_ILA_PPA_ESTIMATIONS_PPA_PROFILE_H__
#define ILANG_ILA_PPA_ESTIMATIONS_PPA_PROFILE_H__

#include <memory>

namespace ilang
{

class PPAProfile_Base
{

public:

    virtual double getBlockTime() = 0;
    virtual double getBlockDynamicPower() = 0;
    virtual double getBlockLeakagePower() = 0;
    virtual double getBlockArea() = 0;
    virtual int getMaximumBitwidth() = 0;

    virtual int getMaximumInstances() = 0;
    virtual bool getIsReusable() = 0;

    // Should not be overriden
    void setGlobalIndex(int index) { m_globalIndex = index; }

    // Should not be overriden
    int getGlobalIndex() { return m_globalIndex;}

private:
    int m_globalIndex;

};

typedef std::shared_ptr<PPAProfile_Base> PPAProfile_ptr;

}

#endif//ILANG_ILA_PPA_ESTIMATIONS_PPA_H__