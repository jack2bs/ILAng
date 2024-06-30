/// \file
/// Implementation of the class PPAAnalyzer

#include "ilang/ila-mngr/u_abs_knob.h"
#include "ilang/ila/ast/expr.h"
#include "ilang/ila/ast/expr_op.h"
#include "ilang/ila/ast_hub.h"
#include <cstddef>
#include <ilang/ppa-estimations/ppa.h>
#include <memory>
#include <unordered_map>

#include "ilang/ppa-estimations/ppa_callbacks.h"
#include "ilang/ppa-estimations/ppa_hardware_block.h"

/// \namespace ilang
namespace ilang {


// This function needs a refactor to remove all of the duplicate code.
// It's super doable just not worth doing yet since it won't actually come
// with a performance improvement, just a sloc one.

double PPAAnalyzer::PerformanceSearchByOperation
(
    const ExprPtr & expr,
    double maximumArgReadyTime,
    PPAAnalysisData & ppaData
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
        
        std::shared_ptr<std::vector<bool>> & memUseTracker =
            ppaData.m_exprToMemUseByCycle.at(expr->arg(0));

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
        double startTime = ((cycleToCheck)  * m_cycleTime);
        return startTime;
    }
    

    HardwareBlock_t hbInd = UidToHardwareBlock(exprOp);
    if (m_blockIsReused[hbInd])
    {
        std::vector<int> & blockTracker = 
            ppaData.m_hardwareBlocksPerCycle[hbInd];

        int maxBlocksPerCycle = (*ppaData.m_maxOfEachBlock)[hbInd];

        // Convert `maximumArgReadyTime` from time to cycles
        size_t firstCycleReady = static_cast<size_t>
            (maximumArgReadyTime / m_cycleTime);

        while (firstCycleReady >= blockTracker.size())
        {
            blockTracker.push_back(0);
        }

        size_t cycleToCheck = firstCycleReady;
        while (blockTracker.at(cycleToCheck) >= maxBlocksPerCycle)
        {
            cycleToCheck++;

            if (cycleToCheck >= blockTracker.size())
            {
                blockTracker.push_back(0);
                break;
            }
        }
        blockTracker.at(cycleToCheck)++;

        double startTime = cycleToCheck * m_cycleTime;
        return startTime;
    }

    return maximumArgReadyTime;

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
        if (ppaData.m_endTimes.find(e) != ppaData.m_endTimes.end())
        {
            return;
        }

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

        ppaData.m_readyTimes.insert({e, maximumArgReadyTime});


        double eProvisionalFinishTime = 0.0;
        double opPerformance = 0.0;
        if (e->is_op())
        {
            SortPtr srt = e->sort();
            int bitwidth = 
                srt->is_bool() ? 1 :
                srt->is_bv() ? srt->bit_width() :
                /*srt->is_mem() ?*/ srt->data_width();

            PPA_Callbacks::s_m_bitwidth = bitwidth;

            // TODO : I don't love the way this deals with memory, think ab it.
            double opPerformance = e->is_mem() ? m_cycleTime :
                m_registrar.getMatchingProfile_LowestBitwidth(
                UidToHardwareBlock(asthub::GetUidExprOp(e)),
                bitwidth)->getBlockTime();

            eProvisionalFinishTime = maximumArgReadyTime + opPerformance;
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

            eProvisionalFinishTime = maximumArgReadyTime;

        }


        /* If it is configured not to pipeline short operations, then we need
         * to wait until the start of the next cycle to start the operation */
        double startCycleUncasted = (maximumArgReadyTime / m_cycleTime);
        int startCycle = static_cast<int>(startCycleUncasted);
        int endCycle = static_cast<int>(eProvisionalFinishTime / m_cycleTime);

        bool spansCycles = startCycle != endCycle;

        bool isShorterThanThreshold = 
            opPerformance < m_configuration.shortOperationLengthThreshold;

        bool shortOpRetime = !m_configuration.pipelineShortOperations 
                                && isShorterThanThreshold;

        bool longOpRetime = m_configuration.pipelineStartAtCycleBoundary;

        double provisionalStartTime = maximumArgReadyTime;
        if ((longOpRetime || shortOpRetime) && spansCycles)
        {
            /* This gets the next cycle except for the case where we are 
             * exactly or almost exactly on a cycle boundary */
            provisionalStartTime = ceil(startCycleUncasted - 0.000001) 
                * m_cycleTime;
        }

        // USEFUL DEBUG DON'T DELETE JACK!
        // ILA_INFO << "Start: " << maximumArgReadyTime << " optime: " << opPerformance << " end: " << eProvisionalFinishTime;

        double startTime = e->is_op()
            ? PerformanceSearchByOperation(e, provisionalStartTime, ppaData) 
            : maximumArgReadyTime;

        ppaData.m_startTimes.insert({e, startTime});

        double eFinishTime = startTime + opPerformance;
        
        ppaData.m_endTimes.insert({e, eFinishTime});

        ILA_ASSERT(eFinishTime - startTime > -0.001) << "finish pre start";

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

    static std::vector<ExprPtr> memVars;

    if (isInit)
    {
        auto end = memVars.end();
        for (const ExprPtr & var : memVars)
        {
            std::shared_ptr<std::vector<bool>> newVec =
                std::make_shared<std::vector<bool>>();

            ppaData.m_exprToMemUseByCycle.insert({var, newVec});
        }

        // TODO: Later, make this based on whether a specific config was 
        // supplied
        ppaData.m_maxOfEachBlock = &m_defaultBlockReuseMax;

        return;
    }

    for (const ExprPtr & var : absknob::GetSttTree(m_ila))
    {
        if (var->is_mem())
        {
            memVars.push_back(var);

            ILA_INFO << "Found mem state, " << var->name().str() 
                << " added vector to list of memVars";
        }
    }

    isInit = true;
    AnalysisDataInitialize(ppaData);
}


/*****************************************************************************/

// This is subject to frequent renaming and changes
void PPAAnalyzer::AnalysisDataDelete(PPAAnalysisData & ppaData)
{
    ppaData.m_exprToMemUseByCycle.clear();
    ppaData.m_endTimes.clear();
    ppaData.m_readyTimes.clear();
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


    int numCycles = static_cast<int>(ceil((ppaData.m_latestTime + 0.00000001)
        / m_cycleTime));



    for (std::vector<int> & vec : ppaData.m_hardwareBlocksPerCycle)
    {
        vec.clear();
        vec.resize(numCycles, 0);
    }
    
    for (auto & pair : ppaData.m_endTimes)
    {

        int count = 0;

        const ExprPtr & expr = pair.first;
        double endTime = pair.second;

        if (!expr->is_op())
        {
            continue;
        }

        AstUidExprOp exprOp = asthub::GetUidExprOp(expr);
        HardwareBlock_t exprHwBlock = UidToHardwareBlock(exprOp);

        double startTime = ppaData.m_startTimes.at(expr);



        int startCycle = static_cast<int>(startTime / m_cycleTime);
        ppaData.m_hardwareBlocksPerCycle[exprHwBlock].at(startCycle)++;
    
    }


    for (size_t i = 0; i < numCycles; i++)
    {
        for (int j = 0; j < m_countBlockTypes; j++)
        {
            std::cout << ' ' << ppaData.m_hardwareBlocksPerCycle[j].at(i);
        }
        std::cout << '\n';
    }
    std::cout << '\n' << '\n';

}

/*****************************************************************************/

const ExprPtr * findDuplicates
(
    const ExprPtr & e,
    std::unordered_map<uint64_t, const ExprPtr> & set
)
{
    if (e->is_var())
    {
        return 0;
    }


    size_t key = 0;
    if (e->is_op())
    {
        for (int i = 0; i < e->arg_num(); i++)
        {
            key |= (e->arg(i)->name().id() & 0x7ffffUL) << (i * 19);
        }
        key |= (asthub::GetUidExprOp(e) & 0x3fUL) << 57UL;
    }
    else if (e->is_const())
    {
        ExprConst & eConst = dynamic_cast<ExprConst&>(*e);
        if (eConst.is_bool())
        {
            key = eConst.val_bool()->val() | (1UL << 62);
        }
        else if (eConst.is_bv())
        {
            key = (eConst.val_bv()->val() << 1) | (1UL << 63);
        }
        else // eConst is mem and we don't expect duplicate constant mems
        {
            return 0;
        }
    }

    while (true)
    {
        auto found = set.find(key);
        if (found == set.end())
        {
            // This path is followed the vast majority of the time
            set.insert({key, e});
            return nullptr;
        }

        const ExprPtr & match = found->second;
        if (match->is_op() && e->is_op())
        {
            /* If the hashes matched and they're both ops, then the ops match
            * as dedicate enough bits to guarantee it. Hence, their arities
            * are equal */
            

            ILA_ASSERT(e->arg_num() == match->arg_num()) << "Arg nums don't match";
            int args = e->arg_num();
            bool theyMatch = true;
            for (int i = 0; i < args; i++)
            {
                theyMatch &=
                    (e->arg(i)->name().id() == match->arg(i)->name().id());
            }
            if (theyMatch)
            {
                return &match;
            }
        
        }
        if (match->is_const() && e->is_const())
        {
            ExprConst & matchConst = dynamic_cast<ExprConst&>(*match);
            ExprConst & eConst = dynamic_cast<ExprConst&>(*e);

            if (matchConst.is_bool() && eConst.is_bool())
            {
                return &match;
            }
            else if (matchConst.is_bv() && eConst.is_bv() 
                && matchConst.val_bv()->val() == eConst.val_bv()->val())
            {
                return &match;
            }
        }
        /* If there was a collision (incorrect match), just increment the
         * key and try again */ 
        static long collisionCount = 0;
        collisionCount++;
        ILA_INFO << "Duplicate search found collision, #" << collisionCount;
        key++;
    }
    // Unreachable
    ILA_ASSERT(false) << "Reached supposedly unreachable code in findDuplicates";
}

/*****************************************************************************/

const ExprPtr * PPAAnalyzer::RemoveDuplicates
(
    const ExprPtr & topExpr,
    std::unordered_map<uint64_t, const ExprPtr> & set
)
{
    for (int i = 0; i < topExpr->arg_num(); i++)
    {
        const ExprPtr & arg = topExpr->arg(i);
        const ExprPtr * newId = RemoveDuplicates(arg, set);
        if (newId)
        {
            // ILA_INFO << "Found duplicate! :)";
            topExpr->replace_arg(i, *newId);
        }
    }
    return findDuplicates(topExpr, set);
}

/*****************************************************************************/

void PPAAnalyzer::PPAAnalyze()
{

    ILA_INFO << "Begin PPA Estimation of " << m_ila;

    std::vector<PPAAnalysisData_ptr> ppaData {};

    for (InstrPtr & instr : absknob::GetInstrTree(m_ila))
    {
        ILA_INFO << "Beginning instruction : " << instr->name().c_str();

        ppaData.emplace_back(std::make_unique<PPAAnalysisData>());
        PPAAnalysisData_ptr & ppaDataTest = ppaData.back();

        AnalysisDataInitialize(*ppaDataTest);

        std::unordered_map<uint64_t, const ExprPtr> set;

        Instr::StateNameSet updated_states = instr->updated_states();
        for (const std::string& s : updated_states) {
            const ExprPtr & update_expr = instr->update(s);

            RemoveDuplicates(update_expr, set);

            PerformanceGet(update_expr, *ppaDataTest, s);
        }
        std::cout << instr->name() << '\n';
        ExtractHardwareBlocks(*ppaDataTest);
    }

    for (PPAAnalysisData_ptr & ppaDataTest : ppaData)
    {
        AnalysisDataDelete(*ppaDataTest);
    }
    

    // ExtractHardwareBlocks(ppaData);
    // AnalysisDataDelete(*ppaDataToTest);
}

/*****************************************************************************/

PPAAnalyzer::PPAAnalyzer(const InstrLvlAbsPtr & ila, double cycle_time)
    : m_ila(ila), m_cycleTime(cycle_time), 
      m_blockIsReused{}, m_defaultBlockReuseMax {}
{

    m_configuration = {
        .pipelineStartAtCycleBoundary = true,

        .pipelineShortOperations = false,
        .shortOperationLengthThreshold = 1.0,

        .reuseShifterBlocks = false,
        .reuseAdditionBlocks = true,
        .reuseMultiplyBlocks = true,
        .reuseDivideBlocks = true,
        .reuseRemainderBlocks = true,

        .maximumShifterBlocks = 2,
        .maximumAdditionBlocks = 2,
        .maximumMultiplyBlocks = 2,
        .maximumDivideBlocks = 2,
        .maximumRemainderBlocks = 2
    };

    ILA_WARN_IF(m_configuration.pipelineStartAtCycleBoundary 
        && m_configuration.pipelineShortOperations) 
        << "Setting pipelineShortOperations while pipelineStartAtCycleBoundary is set will have no effect";

    m_blockIsReused[bShifter] = m_configuration.reuseShifterBlocks;
    m_blockIsReused[bAddition] = m_configuration.reuseAdditionBlocks;
    m_blockIsReused[bMultiplication] = m_configuration.reuseMultiplyBlocks;
    m_blockIsReused[bDivision] = m_configuration.reuseDivideBlocks;
    m_blockIsReused[bRemainder] = m_configuration.reuseRemainderBlocks;
   
    m_defaultBlockReuseMax[bShifter] = m_configuration.maximumShifterBlocks;
    m_defaultBlockReuseMax[bAddition] = m_configuration.maximumAdditionBlocks;
    m_defaultBlockReuseMax[bMultiplication] = m_configuration.maximumMultiplyBlocks;
    m_defaultBlockReuseMax[bDivision] = m_configuration.maximumDivideBlocks;
    m_defaultBlockReuseMax[bRemainder] = m_configuration.maximumRemainderBlocks;   
    
    ILA_ASSERT(m_blockIsReused[bMemory] == false)
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
)   : m_ila(ila), m_cycleTime(cycle_time), m_configuration(*config),
      m_blockIsReused{}, m_defaultBlockReuseMax {}
{

    m_blockIsReused[bShifter] = m_configuration.reuseShifterBlocks;
    m_blockIsReused[bAddition] = m_configuration.reuseAdditionBlocks;
    m_blockIsReused[bMultiplication] = m_configuration.reuseMultiplyBlocks;
    m_blockIsReused[bDivision] = m_configuration.reuseDivideBlocks;
    m_blockIsReused[bRemainder] = m_configuration.reuseRemainderBlocks;
   
    m_defaultBlockReuseMax[bShifter] = m_configuration.maximumShifterBlocks;
    m_defaultBlockReuseMax[bAddition] = m_configuration.maximumAdditionBlocks;
    m_defaultBlockReuseMax[bMultiplication] = m_configuration.maximumMultiplyBlocks;
    m_defaultBlockReuseMax[bDivision] = m_configuration.maximumDivideBlocks;
    m_defaultBlockReuseMax[bRemainder] = m_configuration.maximumRemainderBlocks;  

    ILA_ASSERT(m_blockIsReused[bMemory] == false) 
        << "m_blockReuseConfig not init to 0s";


    ILA_ASSERT(m_configuration.shortOperationLengthThreshold <= cycle_time) 
        << "Short operation threshold is too large";

    registerAllModels(m_registrar);

    PPAAnalyze();
}

}