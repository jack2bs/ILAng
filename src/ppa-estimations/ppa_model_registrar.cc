#include "ilang/ppa-estimations/ppa_hardware_block.h"
#include "ilang/ppa-estimations/ppa_profile_base.h"
#include "ilang/ppa-estimations/ppa_profile_const_example.h"
#include <ilang/ppa-estimations/ppa_model_registrar.h>
#include <memory>
#include <limits.h>
#include <algorithm>

namespace ilang
{

PPA_Registrar::PPA_Registrar()
{
    typedef PPAProfile_Const_Example constProf;
    typedef std::shared_ptr<PPAProfile_Const_Example> constProf_ptr;   
    
    constProf_ptr noHardwareProf {
        std::make_shared<constProf>(0, 0, 0, 0, INT_MAX, -1, false)
    };

    registerProfile(noHardwareProf, bNoHardware);

    /* For now we don't model the memory use power */
    constProf_ptr memoryProf {
        std::make_shared<constProf>(0, 0, 0, 0, INT_MAX, -1, false)
    };

    registerProfile(memoryProf, bMemory);

}

/*****************************************************************************/

void PPA_Registrar::registerProfile
(
    const PPAProfile_ptr & newProfile,
    HardwareBlock_t opType
)
{
    currentlySorted = false;

    m_registeredProfiles.at(opType).push_back(newProfile);

    newProfile->setGlobalIndex(m_numProfiles);
    m_numProfiles++;
}

/*****************************************************************************/

bool sortComparator(const PPAProfile_ptr & profA, const PPAProfile_ptr & profB)
{
    return profA->getMaximumBitwidth() < profB->getMaximumBitwidth();
}

/*****************************************************************************/

void PPA_Registrar::finalizeRegistrar()
{
    for (std::vector<PPAProfile_ptr> & prof : m_registeredProfiles)
    {
        std::sort(prof.begin(), prof.end(), sortComparator);
    }

    currentlySorted = true;
}

/*****************************************************************************/

int PPA_Registrar::getSize()
{
    return m_numProfiles;
}

/*****************************************************************************/

PPAProfile_ptr PPA_Registrar::getMatchingProfile_LowestBitwidth
(
    HardwareBlock_t opType,
    int bitwidth
)
{
    // for lists of the size we expect here, linear search will be faster
    // than binary search
    for (PPAProfile_ptr & prof : m_registeredProfiles.at(opType))
    {
        if (prof->getMaximumBitwidth() >= bitwidth)
        {
            return prof;
        }
    }

    ILA_ERROR << "No profile can handle " << hardwareBlockToString(opType)
        << " with a bitwidth = " << bitwidth;

    // Not reachable without error
    return *(m_registeredProfiles.at(opType).end());
}

/*****************************************************************************/

RegistrarType * PPA_Registrar::getRegisteredProfiles()
{
    return &m_registeredProfiles;
}

}