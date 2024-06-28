/// \file
/// Implementation of the class PPAAnalyzer

#include "ilang/ila-mngr/u_abs_knob.h"
#include "ilang/ila/ast/expr.h"
#include "ilang/ila/ast/expr_op.h"
#include "ilang/ila/ast_hub.h"
#include <ilang/ppa-estimations/ppa.h>
#include <memory>

#include "ilang/ppa-estimations/ppa_callbacks.h"

/// \namespace ilang
namespace ilang {

double PPAAnalyzer::PerformanceSearchByOperation
(
    const ExprPtr & expr,
    double maximumArgReadyTime,
    PPAAnalysisData ppaData
)
{
    /* The methodology below is extra conservative for memory `if_then_else` */
    AstUidExprOp exprOp = asthub::GetUidExprOp(expr);
    if (exprOp == kLoad || exprOp == kStore)
    {
        /* Methodology: Assume for now that each memory device can be written 
         * to or read from within a single cycle, and that only 1 operation is 
         * permited per cycle. Track the memory usage as vectors of bools, 
         * stored in the PPAAnalysisData struct in `m_exprToMemUseByCycle`.
         * Then find and use the earliest available cycle after the operation
         * is ready to execute */

        /* Possible future improvement: Multiple ports on memory; configurable
         * memory access latencies. */

        /* Possible performance issue: Lte and Gte operations aren't primitive,
         * they are implemented instead as (Gt | Eq), and (Lt | Eq) 
         * respectively. This leads to extra operations :( */

        /* Use the mem use tracker associated with the 0th argument which for
         * both loads and stores is the memory state they use */
        
        ILA_INFO << "HERE";
        ILA_INFO << expr->arg(0)->name() << '\t' << expr->arg(0)->is_const()
                 << '\t' << expr->arg(0)->is_mem() << (expr->arg(0)->is_op() ? asthub::GetUidExprOp(expr->arg(0)) : -1);
        std::shared_ptr<std::vector<bool>> memUseTracker = 
            ppaData.m_exprToMemUseByCycle.at(expr->arg(0));

        ILA_INFO << "NOTHERE";

        // Convert `maximumArgReadyTime` from time to cycles
        size_t firstCycleReady = static_cast<size_t>
            (maximumArgReadyTime / m_cycleTime);
        
        while (firstCycleReady >= memUseTracker->size())
        {
            memUseTracker->push_back(false);
        }

        size_t cycleToCheck = firstCycleReady;
        while (memUseTracker->at(cycleToCheck) == true)
        {
            cycleToCheck++;

            if (cycleToCheck >= memUseTracker->size())
            {
                memUseTracker->push_back(true);
                break;
            }
        }

        memUseTracker->at(cycleToCheck) = true;

        // Check this later for off by 1 error (I don't think there's one)
        double finishTime = ((cycleToCheck + 1)  * m_cycleTime);
        return finishTime - maximumArgReadyTime;
    }
    else
    {
        SortPtr srt = expr->sort();
        int bitwidth = 
            srt->is_bool() ? 1 :
            srt->is_bv() ? srt->bit_width() :
            /*srt->is_mem() ?*/ srt->data_width();

        PPA_Callbacks::s_m_bitwidth = bitwidth;
        double time = m_registrar.getMatchingProfile_LowestBitwidth(
            UidToHardwareBlock(exprOp), bitwidth)->getBlockTime();

        return time;
    }
}

/*****************************************************************************/

/* Copied from ilator.cc */
static bool HasLoadFromStore(const ExprPtr& expr) {
    auto monitor = false;
    auto LoadFromStore = [&monitor](const ExprPtr& e) {
        if (e->is_op() && asthub::GetUidExprOp(e) == AstUidExprOp::kLoad) {
            monitor |= e->arg(0)->is_op();
        }
    };
    expr->DepthFirstVisit(LoadFromStore);
    return monitor;
}

/*****************************************************************************/

void PPAAnalyzer::PerformanceGet
(
    const ExprPtr & expr,
    PPAAnalysisData & ppaData,
    const std::string & exprToUpdate
)
{
    if (HasLoadFromStore(expr))
    {
        ILA_ERROR << "There is a load from a store which is not allowed for these estimations";
    }

    /* If `expr` updates a memory state then we need to go through and check 
     * for store operations, to make sure that they're mapped to the vector 
     * which tracks the memory state they store to, which has to be the one 
     * `expr` updates (`exprToUpdate`) */
    if (!exprToUpdate.empty() && expr->is_mem())
    {
        ExprPtr exprToUpdatePtr = nullptr;

        auto findState =
        [&exprToUpdate, &exprToUpdatePtr](const InstrLvlAbsCnstPtr& m) -> void
        {
            ExprPtr temp = m->find_state(exprToUpdate);
            if (temp) {exprToUpdatePtr = temp;}
        };
        m_ila->DepthFirstVisit(findState);

        ILA_INFO << "UPDATING MEM EXPR";

        ILA_ASSERT(ppaData.m_exprToMemUseByCycle.count(exprToUpdatePtr))
            << "Update to memory not in m_exprToMemUseByCycle";

        std::shared_ptr<std::vector<bool>> exprMemoryUsage = 
            ppaData.m_exprToMemUseByCycle.at(exprToUpdatePtr);
        
        std::function<void(const ExprPtr & e)> preCheckForStores =
        [&exprMemoryUsage, &ppaData, &preCheckForStores]
        (const ExprPtr & e)
        {
            if (e->is_mem() && e->is_op())
            {
                ppaData.m_exprToMemUseByCycle.insert(
                    {e, exprMemoryUsage}
                );

                for (size_t i = 0; i < e->arg_num(); i++)
                {
                    preCheckForStores(e->arg(i));
                }
            }
        };

        preCheckForStores(expr);
    }
    auto determineTiming =
    [&ppaData, this](const ExprPtr & e)
    {
        /* If e is a load operation, then we need to make sure that it is
         * mapped to the correct vector for tracking memory usage for the state
         * e loads from */
        if (e->is_op() && asthub::GetUidExprOp(e) == AstUidExprOp::kLoad)
        {
            ILA_ASSERT(ppaData.m_exprToMemUseByCycle.count(e->arg(0)))
                 << "Load op loading from memory not in m_exprToMemUseByCycle";

            std::shared_ptr<std::vector<bool>> exprMemoryUsage = 
                ppaData.m_exprToMemUseByCycle.at(e->arg(0));
                
            ppaData.m_exprToMemUseByCycle.insert(
                {e, exprMemoryUsage}
            );
        }
        /* Methodology: Find time at which all args are ready, and how long the
         * operation will take. It's finish time is the sum of these two times
         */

        double maximumArgReadyTime = 0.0;

        for (size_t i = 0; i < e->arg_num(); i++)
        {
            double currArgReadyTime = ppaData.m_endTimes.at(e->arg(i));
            if (currArgReadyTime > maximumArgReadyTime)
            {
                maximumArgReadyTime = currArgReadyTime; 
            }
        }

        ppaData.m_startTimes.insert({e, maximumArgReadyTime});

        double eReadyTime = 0.0;
        double opPerformance = 0.0;
        if (e->is_op())
        {
            opPerformance = PerformanceSearchByOperation(
                e, maximumArgReadyTime, ppaData
            );
            eReadyTime = maximumArgReadyTime + opPerformance;
        }
        else
        {
            /* If we stumble accross a memory constant, we need to add it to
             * our tracking of memory usage and a list of constant mems so that
             * it can be deleted at the end of the program */
            if (e->is_mem() && e->is_const())
            {
                std::shared_ptr<std::vector<bool>> newVec = 
                     std::make_shared<std::vector<bool>>();
                ppaData.m_exprToMemUseByCycle.insert({e, newVec});
                m_constMems.insert(e);
            }

            // I think this should always be 0 if here
            ILA_WARN_IF(maximumArgReadyTime > 0.001) 
                << "Non op found with start time after 0.0";

            eReadyTime = maximumArgReadyTime;

        }


        /* If it is configured not to pipeline short operations, then we need
         * to wait until the start of the next cycle to start the operation */
        double startCycleUncasted = (maximumArgReadyTime / m_cycleTime);
        int startCycle = static_cast<int>(maximumArgReadyTime / m_cycleTime);
        int endCycle = static_cast<int>(eReadyTime / m_cycleTime);

        bool spansCycles = startCycle != endCycle;

        bool isShorterThanThreshold = 
            opPerformance < m_configuration.shortOperationLengthThreshold;

        bool shortOpRetime = !m_configuration.pipelineShortOperations 
                                && isShorterThanThreshold;

        bool longOpRetime = m_configuration.pipelineStartAtCycleBoundary;

        if ((longOpRetime || shortOpRetime) && spansCycles)
        {
            /* This gets the next cycle except for the case where we are 
             * exactly or almost exactly on a cycle boundary */
            double newStartTime = ceil(startCycleUncasted - 0.000001);
            eReadyTime = ((newStartTime) * m_cycleTime) + opPerformance;
        }

        // USEFUL DEBUG DON'T DELETE JACK!
        // ILA_INFO << "Start: " << maximumArgReadyTime << " optime: " << opPerformance << " end: " << eReadyTime;

        ppaData.m_endTimes.insert({e, eReadyTime});

    };

    expr->DepthFirstVisit(determineTiming);

    double exprEndTime = ppaData.m_endTimes.at(expr);
    if (exprEndTime > ppaData.m_latestTime)
    {
        ppaData.m_latestTime = exprEndTime;
    }

    ILA_INFO << expr->name().str() << " performance estimated at " 
        << ppaData.m_endTimes.at(expr);
}

/*****************************************************************************/

// This is subject to frequent renaming and changes
void PPAAnalyzer::AnalysisDataInitialize(PPAAnalysisData & ppaData)
{
    static bool isInit = false;
    static PPAAnalysisData emptyPpaData = {};

    if (isInit)
    {
        // TODO: Will need to check that this works as intended!
        ppaData = emptyPpaData;
        return;
    }

    for (const ExprPtr & var : absknob::GetSttTree(m_ila))
    {
        if (var->is_mem())
        {
            std::shared_ptr<std::vector<bool>> newVec =
                std::make_shared<std::vector<bool>>();
            ppaData.m_exprToMemUseByCycle.insert({var, newVec});
            emptyPpaData.m_exprToMemUseByCycle.insert({var, newVec});
            ILA_INFO << "Found mem state, " << var->name().str() 
                << " added vector to ppaData at addr " << newVec;
        }
    }

    isInit = true;
}


/*****************************************************************************/

// This is subject to frequent renaming and changes
void PPAAnalyzer::AnalysisDataDelete(PPAAnalysisData & ppaData)
{
    ppaData.m_exprToMemUseByCycle.clear();
    ppaData.m_hardwareBlocksPerCycle.clear();
    ppaData.m_endTimes.clear();
    ppaData.m_startTimes.clear();
    ppaData.m_latestTime = 0.0;
}


/*****************************************************************************/

HardwareBlock_t PPAAnalyzer::UidToHardwareBlock(AstUidExprOp op)
{
    switch (op)
    {
    case kConcatenate:
    case kExtract:
    case kZeroExtend:
    case kSignedExtend:
        return bNoHardware;
    case kNegate:
    case kNot:
    case kComplement:
    case kAnd:
    case kOr:
    case kXor: 
    case kImply:
    case kIfThenElse:
        return bBitwise;
    case kShiftLeft:
    case kArithShiftRight:
    case kLogicShiftRight:
    case kRotateLeft:
    case kRotateRight:
        return bShifter;
    case kAdd:
    case kSubtract:
    case kEqual:
    case kLessThan:
    case kGreaterThan:
    case kUnsignedLessThan:
    case kUnsignedGreaterThan:
        return bAddition;
    case kMultiply:
        return bMultiplication;
    case kDivide:
        return bDivision;
    case kSignedRemainder:
    case kUnsignedRemainder:
    case kSignedModular:
        return bRemainder;
    case kLoad:
    case kStore:
        return bMemory;
    case kApplyFunc:
        return bNoHardware;
    default:
        return bInvalid;
    }
}

/*****************************************************************************/

void PPAAnalyzer::ExtractHardwareBlocks(PPAAnalysisData & ppaData)
{
    /* Methodology: using the ExprPtr's and the ready times in ppaData.
     * m_endTimes, we can determine what hardware blocks are in use during
     * what cycles. */

    int numCycles = static_cast<int>(ceil(ppaData.m_latestTime + 0.00000001
        / m_cycleTime));

    ppaData.m_hardwareBlocksPerCycle.resize(numCycles);
    for (size_t i = 0; i < numCycles; i++)
    {
        for (int j = 0; j < m_countBlockTypes; j++)
        {
            ppaData.m_hardwareBlocksPerCycle.at(i)[j] = 0;
        }
    }
    
    for (auto pair : ppaData.m_endTimes)
    {

        int count = 0;
        for (auto pair2 : ppaData.m_endTimes)
        {
            if (pair2.first == pair.first)
            {
                count++;
            }

            ILA_ASSERT(count <= 1) << "There are duplicates";

        }

        ExprPtr expr = pair.first;
        double endTime = pair.second;

        if (!expr->is_op())
        {
            continue;
        }

        AstUidExprOp exprOp = asthub::GetUidExprOp(expr);
        HardwareBlock_t exprHwBlock = UidToHardwareBlock(exprOp);

        double startTime = ppaData.m_startTimes.at(expr);

        int startCycle = static_cast<int>(startTime / m_cycleTime);

        ppaData.m_hardwareBlocksPerCycle.at(startCycle)[exprHwBlock]++;
    
    }

    // for (size_t i = 0; i < numCycles; i++)
    // {
    //     for (int j = 0; j < m_countBlockTypes; j++)
    //     {
    //         std::cout << ' ' << ppaData.m_hardwareBlocksPerCycle.at(i)[j];
    //     }
    //     std::cout << '\n';
    // }

}

/*****************************************************************************/

void PPAAnalyzer::PPAAnalyze()
{

    ILA_INFO << "Begin PPA Estimation of " << m_ila;

    std::vector<PPAAnalysisData_ptr> ppaData {};


    for (InstrPtr instr : absknob::GetInstrTree(m_ila))
    {
        ILA_INFO << "Beginning instruction : " << instr->name().c_str();

        PPAAnalysisData_ptr ppaDataTest = std::make_shared<PPAAnalysisData>();
        AnalysisDataInitialize(*ppaDataTest);
        ppaData.emplace_back(ppaDataTest);
        
        Instr::StateNameSet updated_states = instr->updated_states();
        for (const std::string& s : updated_states) {
            ExprPtr update_expr = instr->update(s);
            PerformanceGet(update_expr, *ppaDataTest, s);
        }
        ExtractHardwareBlocks(*ppaDataTest);
    }

    for (PPAAnalysisData_ptr ppaDataTest : ppaData)
    {
        AnalysisDataDelete(*ppaDataTest);
    }
    

    // ExtractHardwareBlocks(ppaData);
    // AnalysisDataDelete(*ppaDataToTest);
}

/*****************************************************************************/

PPAAnalyzer::PPAAnalyzer(const InstrLvlAbsPtr & ila, double cycle_time)
    : m_ila(ila), m_cycleTime(cycle_time), m_blockReuseConfig{}
{

    m_configuration = {
        .pipelineStartAtCycleBoundary = true,

        .pipelineShortOperations = false,
        .shortOperationLengthThreshold = 1.0,

        .reuseShifterBlocks = false,
        .reuseAdditionBlocks = false,
        .reuseMultiplyBlocks = true,
        .reuseDivideBlocks = true,
        .reuseRemainderBlocks = true
    };

    ILA_WARN_IF(m_configuration.pipelineStartAtCycleBoundary 
        && m_configuration.pipelineShortOperations) 
        << "Setting pipelineShortOperations while pipelineStartAtCycleBoundary is set will have no effect";

    m_blockReuseConfig[bShifter] = m_configuration.reuseShifterBlocks;
    m_blockReuseConfig[bAddition] = m_configuration.reuseAdditionBlocks;
    m_blockReuseConfig[bMultiplication] = m_configuration.reuseMultiplyBlocks;
    m_blockReuseConfig[bDivision] = m_configuration.reuseDivideBlocks;
    m_blockReuseConfig[bRemainder] = m_configuration.reuseRemainderBlocks;
    ILA_ASSERT(m_blockReuseConfig[bMemory] == false)
        << "m_blockReuseConfig not init to 0s";

    if (m_configuration.shortOperationLengthThreshold > cycle_time)
    {
        m_configuration.shortOperationLengthThreshold = cycle_time;
    }
    

    ILA_ASSERT(m_configuration.shortOperationLengthThreshold <= cycle_time) 
        << "Short operation threshold is too large";

    registerAllModels(m_registrar);

    PPAAnalyze();
}

/*****************************************************************************/

PPAAnalyzer::PPAAnalyzer
(
    const InstrLvlAbsPtr & ila,
    double cycle_time,
    PPAAnalyzerConfig * config
)   :m_ila(ila), m_cycleTime(cycle_time), m_configuration(*config)
{

    m_blockReuseConfig[bShifter] = m_configuration.reuseShifterBlocks;
    m_blockReuseConfig[bAddition] = m_configuration.reuseAdditionBlocks;
    m_blockReuseConfig[bMultiplication] = m_configuration.reuseMultiplyBlocks;
    m_blockReuseConfig[bDivision] = m_configuration.reuseDivideBlocks;
    m_blockReuseConfig[bRemainder] = m_configuration.reuseRemainderBlocks;
    ILA_ASSERT(m_blockReuseConfig[bMemory] == false) 
        << "m_blockReuseConfig not init to 0s";


    ILA_ASSERT(m_configuration.shortOperationLengthThreshold <= cycle_time) 
        << "Short operation threshold is too large";

    registerAllModels(m_registrar);

    PPAAnalyze();
}

}