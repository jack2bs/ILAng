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
#include "ilang/ila/ast/expr_op.h"
#include <cmath>
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

        /* What should the threshold be for a 'short' operation (<= 1.0)
         * default = 1.0 */
        double shortOperationLengthThreshold = 1.0;

        /* The default number of ports on the memory */
        int defaultNumMemoryPorts = 2;


        /* How should implicit register count be estimated? */
        enum regCountMethods 
        {
            /* Slow hardware blocks are adders, shifters, multipliers, 
             * dividers and remainderers*/
            CountSlowHardWareBlocks,

            /* Not recommended can be inaccurate for off critical path times */
            CountCarriedOverCycleBoundaries
        } regCountMethod = CountSlowHardWareBlocks;

        /* Not recommended, significantly overestimates the # of registers */
        bool PutConstantsInRegisters = false;

    };

    PPAAnalyzer
    (
        const InstrLvlAbsPtr & ila,
        double cycle_time,
        const std::string & instr_seq_path,
        const std::string & vcd_path
    );

    PPAAnalyzer
    (
        const InstrLvlAbsPtr & ila,
        double cycle_time,
        const std::string & instr_seq_path,
        const std::string & vcd_path,
        PPAAnalyzerConfig * config
    );

    void PPAAnalyze();

private:

    typedef std::array<std::vector<int>, bNumBlockTypes> UseTracker_t;
    struct PPAAnalysisData
    {
        PPAAnalysisData() : m_latestTime(0) {};

        // The time when each expression is finished in the current analysis
        std::unordered_map<uint64_t, double> m_endTimes;

        // The time when each expression begins in the current analysis
        std::unordered_map<ExprPtr, double> m_startTimes;
        
        // The time when each expression is ready in the current analysis
        std::unordered_map<ExprPtr, double> m_readyTimes;

        // Contains the expression pointers which have been put in registers
        std::unordered_map<ExprPtr, int> m_inRegister;

        // The hardware profile mapped to each expression
        std::unordered_map<ExprPtr, PPAProfile_ptr> m_profiles;

        // Maps each memory expression to it's underlying memory partition's 
        // timing tracker 
        std::unordered_map<ExprPtr, std::shared_ptr<std::vector<int>>> 
            m_exprToMemUseByCycle;

        // Tracks the largest value of any expression in the m_endTimes map
        double m_latestTime;

        // Tracks how many of each block have been used each cycle to enforce
        // configuration maximums
        std::vector<std::vector<int>> m_hardwareUseTracker;

        // Tracker for dynamic programming for the HasLoadFromStore check
        std::unordered_set<const ExprPtr *> m_hasLoadFromStoreVisited;

        //
        // std::array<std::vector<int>, bNumBlockTypes> m_hardwareBlocksPerCycle;
    };
    typedef std::unique_ptr<PPAAnalysisData> PPAAnalysisData_ptr;

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

    /* Returns the time `expr` can begin execution. `maximumArgReadyTime` is 
     * the earliest time at which all args for `expr` can be ready. Uses the 
     * structures in `ppa_time_dat`. */
    double ScheduleAndCount
    (
        const ExprPtr & expr, 
        double maximumArgReadyTime,
        PPAAnalysisData & ppa_time_dat
    );

    /* Initializes `ppa_time_dat` with tracking vectors for all memory states
     * and resizes the hardware block tracking vector with the number of 
     * profiles in the registrar*/
    void AnalysisDataInitialize(PPAAnalysisData & ppa_time_dat);

    /* Cleans up the tracking vectors in `ppa_time_dat` */
    void AnalysisDataDelete(PPAAnalysisData & ppa_time_dat);

    /* Removes duplicates from `topExpr` according to the list in `set`. 
     * `CheckedMap` is used as a visited list (dynamic programming) to speed up
     * the code */
    const ExprPtr * RemoveDuplicates
    (
        const ExprPtr & topExpr, 
        std::unordered_map<uint64_t, const ExprPtr> & set,
        std::unordered_map<uint64_t, const ExprPtr> & checkedMap
    );

    // TODO : downgrading operations with constant operands
    /* Returns the type of hardware block which implements operation `op` */
    HardwareBlock_t UidToHardwareBlock(AstUidExprOp op);

    /* Prints the hardware blocks which are being used (debugging, big 
     * performance hit) */
    void PrintHardwareBlocks(PPAAnalysisData & ppaData);

    /* Fills in m_instrSequence based on data in file : m_instrSeqPath 
     * returns true if m_instrSeqPath could be opened, false otherwise */
    bool MakeInstrSequence();

    /* Count the  # of registers if regCountMethod == CountSlowHardwareBlocks
     * by counting the number of instances of adders, shifters, multipliers,
     * dividers, remainderers, and the number of memory ports */
    void RegCountSlowHwBlocks();

    /* */
    bool MakeVcd();

    /* */
    void FinalizeEstimates
    (
        std::unordered_map<std::string, PPAAnalysisData_ptr> & ppaData,
        bool hadInstrSeq,
        bool hadVcd
    );

    /* */
    void AddToConstDirectory(const ExprPtr & expr);

    const int m_countBlockTypes = HardwareBlock_t::bNumBlockTypes;
    const InstrLvlAbsPtr & m_ila;
    const double m_cycleTime;
    std::set<ExprPtr> m_constMems;
    PPAAnalyzerConfig m_configuration;
    PPA_Registrar m_registrar;
    const std::string & m_instrSeqPath;
    std::vector<std::string> m_instrSequence;
    const std::string & m_vcdPath;
    VCDFile * m_vcdStatistics;
    std::unordered_map<VCDSignalReference, VCDSignalHash*> m_refToHash;
    std::set<uint64_t> m_constDirectory;
    std::vector<ExprPtr> m_memVars;




};

}

#endif//ILANG_ILA_PPA_ESTIMATIONS_PPA_H__