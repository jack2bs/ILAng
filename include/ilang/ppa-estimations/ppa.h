/// \file
/// 

#ifndef ILANG_ILA_PPA_ESTIMATIONS_PPA_H__
#define ILANG_ILA_PPA_ESTIMATIONS_PPA_H__

#include <ilang/ila/instr_lvl_abs.h>
#include <ilang/util/log.h>
#include <ilang/ila/ast/expr.h>
#include <ilang/ila-mngr/u_abs_knob.h>
#include <unordered_map>
#include <vector>
#include <functional>
#include "ilang/ila/ast/expr_op.h"
#include "ilang/ila/ast_hub.h"
#include <cmath>
#include "ilang/ppa-estimations/ppa_profile_base.h"
#include "ppa_hardware_block.h"
#include "ppa_model_registrar.h"
#include <vcdparser/VCDFileParser.hpp>

namespace ilang {

class PPAAnalyzer {

public:

    struct PPAAnalyzerConfig
    {
        /* Should pipelined operations be required to start at a cycle boundary
         * (this includes all operations that get scheduled across cycle 
         * boundaries, which is to say, anything that would need to be 
         * pipelined if not rescheduled, see readme for examples)
         * default = True */
        bool pipelineStartAtCycleBoundary = true;

        /* Should short operations be allowed to span multiple cycles (has no 
         * effect when pipelineStartAtCycleBoundary = true) default = False */
        bool pipelineShortOperations = false;

        /* TODO: Consider changing to be a fraction of cycle_time instead of 
         * absolute time */
        /* What should the threshold be for a 'short' operation (<= cycle_time)
         * default = 1.0 */
        double shortOperationLengthThreshold = 1.0;

        // Which hardware blocks should be reused and which shouldn't
        bool reuseShifterBlocks = false;
        bool reuseAdditionBlocks = false;
        bool reuseMultiplyBlocks = true;
        bool reuseDivideBlocks = true;
        bool reuseRemainderBlocks = true;

        /* How many of each profile of each hardware block can the circuit have
         * (ignored for blocks which aren't marked as reusable, in which case 
         * the required number of the block is instantiated) */

        // TODO: Implement the use of these configuration parameters
        int maximumShifterBlocks;
        int maximumAdditionBlocks;
        int maximumMultiplyBlocks;
        int maximumDivideBlocks;
        int maximumRemainderBlocks;

    };

    PPAAnalyzer(const InstrLvlAbsPtr & ila, double cycle_time);

    PPAAnalyzer
    (
        const InstrLvlAbsPtr & ila,
        double cycle_time,
        PPAAnalyzerConfig * config
    );

    void PPAAnalyze();

private:

    typedef std::array<std::vector<int>, bNumBlockTypes> UseTracker_t;
    struct PPAAnalysisData
    {
        PPAAnalysisData() : m_latestTime(0) {};

        // The time when each expression is finished in the current analysis
        std::unordered_map<ExprPtr, double> m_endTimes;

        // The time when each expression begins in the current analysis
        std::unordered_map<ExprPtr, double> m_startTimes;
        
        // The time when each expression is ready in the current analysis
        std::unordered_map<ExprPtr, double> m_readyTimes;

        // Maps each memory expression to it's underlying memory partition's 
        // timing tracker 
        std::unordered_map<ExprPtr, std::shared_ptr<std::vector<bool>>> 
            m_exprToMemUseByCycle;

        // Tracks the largest value of any expression in the m_endTimes map
        double m_latestTime;

        // Tracks how many of each block have been used each cycle to enforce
        // configuration maximums
        std::vector<std::vector<int>> m_hardwareUseTracker;

        //
        // std::array<std::vector<int>, bNumBlockTypes> m_hardwareBlocksPerCycle;
    };
    typedef std::unique_ptr<PPAAnalysisData> PPAAnalysisData_ptr;

    // struct HardwareBlockCounter
    // {
    //     int bitwiseOps;
    //     int additionOps;
    //     int multiplyOps;
    //     int divideOps;
    //     int remainderOps;
    // };

    /* Estimates the performance of `expr`, which updates the optional 
     * parameter `exprToUpdate`, and uses the structures in `ppa_time_dat`, 
     * which should be reset outside of this function when applicable 
     * (switching between non-simultaneous expressions) */
    void PerformanceGet
    (
        const ExprPtr & expr,
        PPAAnalysisData & ppa_time_dat,
        const std::string & exprToUpdate = ""
    );

    /* Returns the time between `expr`'s arguments being ready and it's 
     * execution being finished. `maximumArgReadyTime` is the earliest time at 
     * which `expr` can start, and uses the structures in `ppa_time_dat`. Fills
     * in start times */
    double ScheduleAndCount
    (
        const ExprPtr & expr, 
        double maximumArgReadyTime,
        PPAAnalysisData & ppa_time_dat
    );

    /* Initializes `ppa_time_dat` with tracking vectors for all memory states
     */
    void AnalysisDataInitialize(PPAAnalysisData & ppa_time_dat);

    /* Cleans up the tracking vectors in `ppa_time_dat` */
    void AnalysisDataDelete(PPAAnalysisData & ppa_time_dat);

    /* Removes duplicates from the expressions according to the list */
    const ExprPtr * RemoveDuplicates
    (
        const ExprPtr & topExpr, 
        std::unordered_map<uint64_t, const ExprPtr> & set
    );


    /* Returns the type of hardware block which implements operation `op` */
    HardwareBlock_t UidToHardwareBlock(AstUidExprOp op);
    /* Prints the hardware blocks which are being used */
    void PrintHardwareBlocks(PPAAnalysisData & ppaData);

    const int m_countBlockTypes = HardwareBlock_t::bNumBlockTypes;
    const InstrLvlAbsPtr & m_ila;
    const double m_cycleTime;
    std::set<ExprPtr> m_constMems;
    PPAAnalyzerConfig m_configuration;
    std::array<bool,bNumBlockTypes> m_blockIsReused;
    std::array<int,bNumBlockTypes> m_defaultBlockReuseMax;
    PPA_Registrar m_registrar;
    VCDFile * m_vcdStatistics;

};

}

#endif//ILANG_ILA_PPA_ESTIMATIONS_PPA_H__