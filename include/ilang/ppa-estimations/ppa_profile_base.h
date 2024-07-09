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

    // Should not be overriden
    void incNumInstances() { m_numInstances++; }

    // Should not be overriden
    void setNumInstances(int newNum) { m_numInstances = newNum; }

    // Should not be overriden
    int getNumInstances() { return m_numInstances; }

private:
    int m_globalIndex = 0;
    int m_numInstances = 0;

};

typedef std::shared_ptr<PPAProfile_Base> PPAProfile_ptr;

}

#endif//ILANG_ILA_PPA_ESTIMATIONS_PPA_H__