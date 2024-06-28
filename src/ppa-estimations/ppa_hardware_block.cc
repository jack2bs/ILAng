#include <ilang/ppa-estimations/ppa_hardware_block.h>

namespace ilang
{

std::string hardwareBlockToString(HardwareBlock_t op)
{
    std::string names [HardwareBlock_t::bNumBlockTypes] {
        "no hardware",
        "bitwise operation",
        "+ (adder type)",
        "* (multiplier type)",
        "/ (division type)",
        "Memory operation",
        "Remainder operation",
    };
    return names[op];
}

}