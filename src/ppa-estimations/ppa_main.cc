/// \file
/// Implementation of the class PPAAnalyzer

#include "VCDParser.hpp"
#include "ilang/ila-mngr/u_abs_knob.h"
#include "ilang/ila/ast/expr.h"
#include "ilang/ila/ast/expr_op.h"
#include "ilang/ila/ast_hub.h"
#include <chrono>
#include <cstddef>
#include <fstream>
#include <ilang/ppa-estimations/ppa.h>
#include <memory>
#include <unordered_map>

#include "ilang/ppa-estimations/ppa_callbacks.h"
#include "ilang/ppa-estimations/ppa_hardware_block.h"
#include "ilang/ppa-estimations/ppa_profile_base.h"

/// \namespace ilang
namespace ilang {


void PPAAnalyzer::AddToConstDirectory(const ExprPtr & expr)
{
    if (!m_configuration.PutConstantsInRegisters)
    {
        return;
    }

    uint64_t id = expr->name().id();
    if (m_constDirectory.find(id) != m_constDirectory.end())
    {
        return;
    }
    
    m_constDirectory.insert(id);

    m_registrar.getMatchingProfile_LowestBitwidth(
        bRegisterHold, 
        (expr->is_bool() ? 1 : expr->sort()->bit_width())
    )->incNumInstances();
}

/*****************************************************************************/

// This function needs a refactor to remove all of the duplicate code.
// It's super doable just not worth doing yet since it won't actually come
// with a performance improvement, just a sloc one.

double PPAAnalyzer::ScheduleAndCount
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
         * to or read from in one whole cycle, and that only 1 operation is 
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

        std::shared_ptr<std::vector<int>> & memUseTracker =
            ppaData.m_exprToMemUseByCycle.at(expr->arg(0));

        // Convert `maximumArgReadyTime` from time to cycles
        size_t firstCycleReady = static_cast<size_t>
            (maximumArgReadyTime / m_cycleTime);
        
        while (firstCycleReady >= memUseTracker->size())
        {
            memUseTracker->push_back(0);
        }

        size_t cycleToCheck = firstCycleReady;
        while (memUseTracker->at(cycleToCheck) >= 
            m_configuration.defaultNumMemoryPorts)
        {
            cycleToCheck++;

            if (cycleToCheck >= memUseTracker->size())
            {
                memUseTracker->push_back(0);
                break;
            }
        }

        memUseTracker->at(cycleToCheck)++;

        // Check this later for off by 1 error (I don't think there's one)
        double startTime = ((cycleToCheck)  * m_cycleTime);
        return startTime;
    }


    HardwareBlock_t hbInd = UidToHardwareBlock(exprOp);

    const SortPtr & srt = expr->sort();
    int bitwidth = 
        srt->is_bool() ? 1 :
        srt->is_bv() ? srt->bit_width() :
        /*srt->is_mem() ?*/ srt->data_width();

    const PPAProfile_ptr & prof = 
        m_registrar.getMatchingProfile_LowestBitwidth(hbInd, bitwidth);

    int profileIndex = prof->getGlobalIndex();

    // std::vector<int> & blockTracker =
    //     ppaData.m_hardwareBlocksPerCycle[profileIndex];

    std::vector<int> & blockTracker = 
        ppaData.m_hardwareUseTracker.at(profileIndex);

    int maxBlocksPerCycle = prof->getMaximumInstances();

    // Convert `maximumArgReadyTime` from time to cycles
    size_t firstCycleReady = static_cast<size_t>
        (maximumArgReadyTime / m_cycleTime);

    while (firstCycleReady >= blockTracker.size())
    {
        blockTracker.push_back(0);
    }

    size_t cycleToCheck = firstCycleReady;

    double startTime = maximumArgReadyTime;

    if (prof->getIsReusable())
    {
        while (blockTracker.at(cycleToCheck) >= maxBlocksPerCycle)
        {
            cycleToCheck++;

            if (cycleToCheck >= blockTracker.size())
            {
                blockTracker.push_back(0);
                break;
            }
        }

        if (prof->getNumInstances() == blockTracker.at(cycleToCheck))
        {
            prof->incNumInstances();
        }
    }
    else
    {
        prof->incNumInstances();
    }

    blockTracker.at(cycleToCheck)++;
    startTime = fmax(cycleToCheck * m_cycleTime, startTime);

    return startTime;
}

/*****************************************************************************/

/* Adapted from ilator.cc - added dynamic programming :) */
static bool HasLoadFromStore
(
    const ExprPtr& expr,
    std::unordered_set<const ExprPtr *> & visited
) 
{
    auto monitor = false;

    std::function<void(const ExprPtr & e)> LoadFromStore = 
    [&monitor, &visited, &LoadFromStore] (const ExprPtr& e) {
        if (visited.count(&e))
        {
            return;
        }

        if (e->is_op() && asthub::GetUidExprOp(e) == AstUidExprOp::kLoad) {
            monitor |= e->arg(0)->is_op();
        }
        visited.insert(&e);

        for (int i = 0; i < e->arg_num(); i++)
        {
            LoadFromStore(e->arg(i));
        }
    };
    LoadFromStore(expr);
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

    if (HasLoadFromStore(expr, ppaData.m_hasLoadFromStoreVisited))
    {
        ILA_ERROR << "There is a load from a store which is not allowed for these estimations";
    }

    /* If there is a constant, and it's not a memory constant, then add it to 
     * the const directory. This const directory is currently used in the case 
     * that the configuration wants constants to be put in registers. The
     * directory makes sure that duplicate registers aren't made for duplicate
     * constants. This configuration setting is highly discouraged as constants
     * are often used to represent ideas that wouldn't really be programmed as 
     * constants. */
    if (expr->is_const() && !expr->is_mem())
    {
        AddToConstDirectory(expr);
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

        std::shared_ptr<std::vector<int>> exprMemoryUsage = 
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

                if (asthub::GetUidExprOp(e) == kStore)
                {
                    preCheckForStores(e->arg(0));
                }
                else if (asthub::GetUidExprOp(e))
                {
                    preCheckForStores(e->arg(1));
                    preCheckForStores(e->arg(2));
                }
                

                // for (size_t i = 0; i < e->arg_num(); i++)
                // {

                //     preCheckForStores(e->arg(i));
                // }
            }
        };

        preCheckForStores(expr);
    }

    auto determineTiming =
    [&ppaData, this](const ExprPtr & e)
    {

        /* If e is a constant, then we should add it to the const directory
         * if it's not there already. The const directory is used for counting
         * the number of consts in case they need to be put in registers. */


        /* If e is a load operation, then we need to make sure that it is
         * mapped to the correct vector for tracking memory usage for the state
         * e loads from */
        if (e->is_op() && asthub::GetUidExprOp(e) == AstUidExprOp::kLoad)
        {
            ILA_ASSERT(ppaData.m_exprToMemUseByCycle.count(e->arg(0)))
                 << "Load op loading from memory not in m_exprToMemUseByCycle";

            std::shared_ptr<std::vector<int>> exprMemoryUsage = 
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
            double currArgReadyTime = ppaData.m_endTimes.at(e->arg(i)->name().id());
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

            AstUidExprOp exprOp = asthub::GetUidExprOp(e);

            PPAProfile_ptr prof = m_registrar.getMatchingProfile_LowestBitwidth
                (UidToHardwareBlock(exprOp), bitwidth);

            // TODO : I don't love the way this deals with memory, think ab it.

            if (exprOp == kLoad)
            {
                opPerformance = m_cycleTime;
            }
            else if (exprOp == kStore)
            {
                opPerformance = 0.0;
            }
            else
            {
                opPerformance = prof->getBlockTime();
            }

            ppaData.m_profiles.insert({e, prof});

            eProvisionalFinishTime = maximumArgReadyTime + opPerformance;
        }
        else
        {
            /* If we stumble accross a memory constant, we need to add it to
             * our tracking of memory usage and a list of constant mems so that
             * it can be deleted at the end of the program */
            if (e->is_mem() && e->is_const())
            {
                std::shared_ptr<std::vector<int>> newVec = 
                     std::make_shared<std::vector<int>>();
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
        double readyCycleUncasted = (maximumArgReadyTime / m_cycleTime);
        int readyCycle = static_cast<int>(readyCycleUncasted);
        int endCycle = static_cast<int>(eProvisionalFinishTime / m_cycleTime);

        bool spansCycles = readyCycle != endCycle;

        bool isShorterThanThreshold = opPerformance 
            < (m_configuration.shortOperationLengthThreshold * m_cycleTime);

        bool shortOpRetime = !m_configuration.pipelineShortOperations 
                                && isShorterThanThreshold;

        bool longOpRetime = m_configuration.pipelineStartAtCycleBoundary;

        double provisionalStartTime = maximumArgReadyTime;
        if ((longOpRetime || shortOpRetime) && spansCycles)
        {
            /* This gets the next cycle except for the case where we are 
             * exactly or almost exactly on a cycle boundary */
            provisionalStartTime = ceil(readyCycleUncasted - 0.000001)
                * m_cycleTime;
        }

        double startTime = e->is_op()
            ? ScheduleAndCount(e, provisionalStartTime, ppaData) 
            : maximumArgReadyTime;

        ppaData.m_startTimes.insert({e, startTime});

        size_t startCycle = static_cast<size_t>(startTime / m_cycleTime);

        /* Determine register usage */
        int num_args = e->arg_num();
        for (int i = 0; i < num_args; i++)
        {
            if (e->arg(i)->is_mem() || e->arg(i)->is_const())
            {
                continue;
            }

            SortPtr srt = e->arg(i)->sort();
            int bitwidth = 
                srt->is_bool() ? 1 :
                srt->is_bv() ? srt->bit_width() :
                /*srt->is_mem() ?*/ srt->data_width();

            const PPAProfile_ptr & holdProf = m_registrar.
                getMatchingProfile_LowestBitwidth(bRegisterHold, bitwidth);
            int holdProfInd = holdProf->getGlobalIndex();
            
            const PPAProfile_ptr & readProf = m_registrar.
                getMatchingProfile_LowestBitwidth(bRegisterRead, bitwidth);
            int readProfInd = readProf->getGlobalIndex();

            const PPAProfile_ptr & writeProf = m_registrar.
                getMatchingProfile_LowestBitwidth(bRegisterWrite, bitwidth);
            int writeProfInd = writeProf->getGlobalIndex();


            int argReadyCycle = static_cast<size_t>
                (ppaData.m_endTimes.at(e->arg(i)->name().id()) / m_cycleTime);

            /* Extract registers - if data spans cycles then it needs to 
             * be put in a register. Writing to and reading from a register, 
             * for now, does not come with a delay. Delays specified in the 
             * profiles will be ignored. They do consume power and area. 
             * leakage power and area are specified in the hold profile, and 
             * dynamic power is specified in the read and write profiles. A 
             * maximum number of registers may not be specified, as the number
             * must be implied.
             * 
             * This definitely makes sense, in my opinion, for counting reg 
             * reads and writes. I think in terms of countning the total number
             * of registers, this heurisitic should just change to the number 
             * slow hardware blocks. Possibly #add + #shift + #mul etc. TODO:
             * consider this decision. */
            if (argReadyCycle != startCycle)
            {
                if (ppaData.m_inRegister.find(e->arg(i)) 
                    == ppaData.m_inRegister.end())
                {
                    std::vector<int> & writeVec = 
                        ppaData.m_hardwareUseTracker.at(writeProfInd);

                    if (writeVec.size() <= argReadyCycle)
                    {
                        writeVec.resize(argReadyCycle+1, 0);
                    }

                    int writes = writeVec.at(argReadyCycle)++;

                    if (!writeProf->getIsReusable() 
                        || writeProf->getNumInstances() < writes)
                    {
                        writeProf->incNumInstances();
                    }

                    ppaData.m_inRegister.insert({e->arg(i), 0});
                }

                if ((ppaData.m_inRegister.at(e->arg(i)) < startCycle )
                 && (m_configuration.regCountMethod 
                 == PPAAnalyzerConfig::CountCarriedOverCycleBoundaries))
                {
                    std::cout << "HERE\n";

                    int heldTill = std::max(
                        argReadyCycle, ppaData.m_inRegister.at(e->arg(i))
                    );

                    std::vector<int> & holdVec =
                        ppaData.m_hardwareUseTracker.at(holdProfInd);

                    if (holdVec.size() <= startCycle)
                    {
                        holdVec.resize(startCycle+1, 0);
                    }

                    for (; heldTill < startCycle; heldTill++)
                    {
                        int holds = holdVec.at(heldTill)++;
                        
                        if (!holdProf->getIsReusable()
                            || holdProf->getNumInstances() < holds)
                        {
                            holdProf->incNumInstances();
                        }
                    }

                    ppaData.m_inRegister.at(e->arg(i)) = startCycle;
                }

                std::vector<int> & readVec =
                    ppaData.m_hardwareUseTracker.at(readProfInd);

                if (readVec.size() <= startCycle)
                {
                    readVec.resize(startCycle+1, 0);
                }

                int reads = readVec.at(startCycle)++;
                
                if (!readProf->getIsReusable()
                    || readProf->getNumInstances() < reads)
                {
                    readProf->incNumInstances();
                }
            }

        }


        double eFinishTime = startTime + opPerformance;
        
        ppaData.m_endTimes.insert({e->name().id(), eFinishTime});

        ILA_ASSERT(eFinishTime - startTime > -0.001) << "finish pre start";

    };

    auto preVisit =
    [&ppaData, this](const ExprPtr & e)
    {

        if (ppaData.m_endTimes.find(e->name().id()) != ppaData.m_endTimes.end())
        {
            return true;
        }
        return false;
    };

    std::function<void(const ExprPtr & e)> dfvpp =
    [&ppaData, &preVisit, &determineTiming, &dfvpp](const ExprPtr & e) -> void
    {

        if (ppaData.m_endTimes.find(e->name().id()) != ppaData.m_endTimes.end())
        {
            return;
        }
        for (int i = 0; i < e->arg_num(); i++)
        {
            dfvpp(e->arg(i));
        }
        determineTiming(e);
    };

    dfvpp(expr);


    // expr->DepthFirstVisit(determineTiming);

    double exprEndTime = ppaData.m_endTimes.at(expr->name().id());
    if (exprEndTime > ppaData.m_latestTime)
    {
        ppaData.m_latestTime = exprEndTime;
    }

    ILA_INFO << expr->name().str() << " performance estimated at " 
        << ppaData.m_endTimes.at(expr->name().id());
}

/*****************************************************************************/

// This is subject to frequent renaming and changes
void PPAAnalyzer::AnalysisDataInitialize(PPAAnalysisData & ppaData)
{
    static bool isInit = false;
    static PPAAnalysisData emptyPpaData = {};

    if (isInit)
    {
        auto end = m_memVars.end();
        for (const ExprPtr & var : m_memVars)
        {
            std::shared_ptr<std::vector<int>> newVec =
                std::make_shared<std::vector<int>>();

            ppaData.m_exprToMemUseByCycle.insert({var, newVec});
        }

        ppaData.m_hardwareUseTracker.resize(m_registrar.getSize());

        return;
    }

    for (const ExprPtr & var : absknob::GetSttTree(m_ila))
    {
        if (var->is_mem())
        {
            m_memVars.push_back(var);

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
    case kEqual:
        return bBitwise;
    case kShiftLeft:
    case kArithShiftRight:
    case kLogicShiftRight:
    case kRotateLeft:
    case kRotateRight:
        return bShifter;
    case kAdd:
    case kSubtract:
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

void PPAAnalyzer::PrintHardwareBlocks(PPAAnalysisData & ppaData)
{

    size_t numCycles = static_cast<size_t>(
        ceil((ppaData.m_latestTime + 0.0000001) / m_cycleTime)
    );

    size_t registrarSize = m_registrar.getSize();

    for (size_t i = 0; i < numCycles; i++)
    {
        for (int j = 0; j < registrarSize; j++)
        {
            int numToPrint = (ppaData.m_hardwareUseTracker.at(j).size() <= i) 
                ? 0 : ppaData.m_hardwareUseTracker.at(j).at(i);

            std::cout << ' ' << numToPrint;
        }

        for (auto state : absknob::GetSttTree(m_ila))
        {
            if (state->is_mem())
            {
                int numToPrint = (ppaData.m_exprToMemUseByCycle.at(state)->size() <= i)
                    ? 0 : ppaData.m_exprToMemUseByCycle.at(state)->at(i);

                std::cout << ' ' << numToPrint;
            }
            
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
    std::unordered_map<uint64_t, const ExprPtr> & set,
    std::unordered_map<uint64_t, const ExprPtr> & checkedMap
)
{
    for (int i = 0; i < topExpr->arg_num(); i++)
    {
        const ExprPtr & arg = topExpr->arg(i);
        const ExprPtr * newId = nullptr;

        auto found = checkedMap.find(arg->name().id());
        if (found != checkedMap.end())
        {
            newId = &found->second;
        }
        else
        {
            const ExprPtr * newId = RemoveDuplicates(arg, set, checkedMap);
            if (newId)
            {
                checkedMap.insert({arg->name().id(), *newId});
            }
        }


        if (newId)
        {
            // ILA_INFO << "Found duplicate! :)";
            topExpr->replace_arg(i, *newId);
        }
    }
    return findDuplicates(topExpr, set);
}

/*****************************************************************************/

bool PPAAnalyzer::MakeInstrSequence()
{
    if (!m_instrSeqPath.compare(""))
    {
        return false;
    }

    std::ifstream instr_seq_in(m_instrSeqPath, std::ios_base::in);

    if (!instr_seq_in.is_open())
    {
        ILA_WARN << "Name of instruction sequence file: " << m_instrSeqPath <<
            "could not be found/opened";

        return false;
    }
    
    std::string instrLine;
    size_t instrCount;

    while (std::getline(instr_seq_in, instrLine))
    {
        std::string & instrName = m_instrSequence.emplace_back();

        instrName = instrLine.substr(instrLine.find('\t') + 1);
        instrName = instrName.substr(0, instrName.find(' '));

        instrCount++;
    }

    return true;
}

/*****************************************************************************/

bool PPAAnalyzer::MakeVcd()
{
    VCDFileParser parser = VCDFileParser();
    m_vcdStatistics = parser.parse_file(m_vcdPath);
    if (!m_vcdStatistics)
    {
        return false;
    }
    

    std::vector<VCDSignal *> * signals = m_vcdStatistics->get_signals();
    for (VCDSignal * sig : *signals)
    {
        m_refToHash.insert({sig->reference, &sig->hash});
    }

    return true;
}


/*****************************************************************************/

void PPAAnalyzer::RegCountSlowHwBlocks()
{
    int numRegs = m_memVars.size();
    for (const ExprPtr & memVar : m_memVars)
    {
        int bw = memVar->sort()->data_width();
        PPAProfile_ptr regProf = m_registrar
            .getMatchingProfile_LowestBitwidth(bRegisterHold, bw);
        
        regProf->setNumInstances(
            regProf->getNumInstances() + m_configuration.defaultNumMemoryPorts
        );
    }

    RegistrarType & regis = *m_registrar.getRegisteredProfiles();

    int hwb = 0;
    for (const std::vector<PPAProfile_ptr> & vec : regis)
    {
        static std::array<bool, bNumBlockTypes> slowHardwareBlock = 
            {false, false, true, true, true, true, 
             false, true, false, false, false };

        if (slowHardwareBlock[hwb])
        {
            for (const PPAProfile_ptr & prof : vec)
            {
                int bw = prof->getMaximumBitwidth();
                PPAProfile_ptr regProf = m_registrar
                    .getMatchingProfile_LowestBitwidth(bRegisterHold, bw);

                std::cout << regProf->getNumInstances() << ' ' << prof->getNumInstances() << '\n';
                
                regProf->setNumInstances( 
                    regProf->getNumInstances() + prof->getNumInstances()
                );

                numRegs += prof->getNumInstances();
            }
        }
        hwb++;
    }
}

/*****************************************************************************/

void PPAAnalyzer::FinalizeEstimates
(
    std::unordered_map<std::string, PPAAnalysisData_ptr> & ppaData,
    bool hadInstrSeq,
    bool hadVcd
)
{
    std::ofstream outFile("PPAOutput", std::ios_base::trunc);


    /* Count the  # of registers if regCountMethod == CountSlowHardwareBlocks
     * by counting the number of instances of adders, shifters, multipliers,
     * dividers, remainderers, and the number of memory ports */
    if (m_configuration.regCountMethod == 
        PPAAnalyzerConfig::CountSlowHardWareBlocks)
    {
        RegCountSlowHwBlocks();
    }

    RegistrarType & regis = *m_registrar.getRegisteredProfiles();

    outFile << "Leakage Power Information:" << std::endl << std::endl;

    // Start by getting hardware block counts for leakage and area
    double runningLeakageCount = 0.0;
    for (const std::vector<PPAProfile_ptr> & vec : regis)
    {
        for (const PPAProfile_ptr & prof : vec)
        {
            double newLeakage = 
                prof->getNumInstances() * prof->getBlockLeakagePower();

            outFile << newLeakage << std::endl;
            runningLeakageCount += newLeakage;
        }
    }

    outFile << "Total Leakage power: " << runningLeakageCount 
        << std::endl << std::endl;

    outFile << "Area Information:" << std::endl << std::endl;

    double runningAreaCount = 0.0;
    for (const std::vector<PPAProfile_ptr> & vec : regis)
    {
        for (const PPAProfile_ptr & prof : vec)
        {
            double newArea = 
                prof->getNumInstances() * prof->getBlockArea();

            outFile << newArea << std::endl;
            runningAreaCount += newArea;
        }
    }
    outFile << "Total area: " << runningAreaCount << std::endl << std::endl;


    if (!hadInstrSeq && !hadVcd)
    {

        outFile << "Timing Information:" << std::endl << std::endl;
        for (auto & data : ppaData)
        {
            PPAAnalysisData_ptr & ppaDataAnalysis = data.second;

            outFile << data.first << " : " << ppaDataAnalysis->m_latestTime;
            outFile << " (" << static_cast<int>(
                ceil(ppaDataAnalysis->m_latestTime / m_cycleTime));
            outFile << "cycles)";
            outFile << std::endl;
        }

        outFile << std::endl;

        PPA_Callbacks::s_m_currSwitchingActivity = 0.5;

        outFile << "Dynamic Power Information: " << std::endl << std::endl;
        for (auto & data : ppaData)
        {
            PPAAnalysisData_ptr & ppaDataAnalysis = data.second;

            double instrDynPwr = 0.0;
            for (auto & pair : ppaDataAnalysis->m_startTimes)
            {
                const ExprPtr & expr = pair.first;
                if (!expr->is_op())
                {
                    continue;
                }


                const SortPtr & srt = expr->sort();
                int bitwidth = 
                    srt->is_bool() ? 1 :
                    srt->is_bv() ? srt->bit_width() :
                    /*srt->is_mem() ?*/ srt->data_width();

                PPA_Callbacks::s_m_bitwidth = bitwidth;
                
                instrDynPwr += ppaDataAnalysis->m_profiles.at(expr)
                    ->getBlockDynamicPower();
            }

            double numCycles = 
                ceil((ppaDataAnalysis->m_latestTime + 0.000001) / m_cycleTime);


            outFile << data.first << " : " << instrDynPwr / numCycles;
            outFile << std::endl;
        }
    }
    else if (hadInstrSeq && !hadVcd)
    {
        outFile << "TODO";
    }
    else if (hadInstrSeq && hadVcd)
    {
        outFile << "TODO";
    }

    ILA_INFO << "DONE!!!";

    outFile.close();
    
}

/*****************************************************************************/

void PPAAnalyzer::PPAAnalyze()
{

    ILA_INFO << "Begin PPA Estimation of " << m_ila;

    // std::vector<PPAAnalysisData_ptr> ppaData {};

    std::unordered_map<std::string, PPAAnalysisData_ptr> ppaData {};

    std::unordered_map<uint64_t, const ExprPtr> set;
    std::unordered_map<uint64_t, const ExprPtr> checkedMap;

    {// for repeated naming purposes

    ILA_INFO << "Beginning fetch, valid, decode";

    // Wow I hate this syntax.
    PPAAnalysisData_ptr & ppaDataTest = 
        ppaData.insert(
            {"__fvd__", std::make_unique<PPAAnalysisData>()}
        ).first->second;

    AnalysisDataInitialize(*ppaDataTest);



    auto PerIla = 
    [this, &set, &checkedMap, &ppaDataTest](const InstrLvlAbsCnstPtr & m)
    {
        const ExprPtr & fetch_expr = m->fetch();
        if (fetch_expr)
        {        
            RemoveDuplicates(fetch_expr, set, checkedMap);
            PerformanceGet(fetch_expr, *ppaDataTest);
        }

        const ExprPtr & valid_expr = m->valid();
        if (valid_expr)
        {
            RemoveDuplicates(valid_expr, set, checkedMap);
            PerformanceGet(valid_expr, *ppaDataTest);
        }
    };
    m_ila->DepthFirstVisit(PerIla);

    for (InstrPtr & instr : absknob::GetInstrTree(m_ila))
    {
        const ExprPtr & decode_expr = instr->decode();
        RemoveDuplicates(decode_expr, set, checkedMap);
        PerformanceGet(decode_expr, *ppaDataTest);
    }

    ppaDataTest->m_hasLoadFromStoreVisited.clear();

    }

    for (InstrPtr & instr : absknob::GetInstrTree(m_ila))
    {
        ILA_INFO << "Beginning instruction : " << instr->name().c_str();

        PPAAnalysisData_ptr & ppaDataTest = ppaData.insert(
            {instr->name().str(), std::make_unique<PPAAnalysisData>()}
        ).first->second;

        AnalysisDataInitialize(*ppaDataTest);

        // std::unordered_map<uint64_t, const ExprPtr> set;
        // std::unordered_map<uint64_t, const ExprPtr> checkedMap;

        Instr::StateNameSet updated_states = instr->updated_states();
        for (const std::string& s : updated_states) {
            const ExprPtr & update_expr = instr->update(s);

            RemoveDuplicates(update_expr, set, checkedMap);
            PerformanceGet(update_expr, *ppaDataTest, s);
        }

        // PrintHardwareBlocks(*ppaDataTest);

        // ExtractHardwareBlocks(*ppaDataTest);


        ppaDataTest->m_hasLoadFromStoreVisited.clear();
    }

    // for (auto & ppaDataTest : ppaData)
    // {
    //     AnalysisDataDelete(*ppaDataTest.second);
    // }

    bool hadInstrSeq = MakeInstrSequence();
    bool hadVcd = false;
    if (hadInstrSeq)
    {
        hadVcd = MakeVcd();
    }

    FinalizeEstimates(ppaData, hadInstrSeq, hadVcd);

    // ExtractHardwareBlocks(ppaData);
    // AnalysisDataDelete(*ppaDataToTest);
}

/*****************************************************************************/

PPAAnalyzer::PPAAnalyzer
(
    const InstrLvlAbsPtr & ila,
    double cycle_time,
    const std::string & instr_seq_path,
    const std::string & vcd_path
)
    : m_ila(ila), m_cycleTime(cycle_time), m_instrSeqPath(instr_seq_path),
      m_vcdPath(vcd_path)
{

    m_configuration = {
        .pipelineStartAtCycleBoundary = true,

        .pipelineShortOperations = false,
        .shortOperationLengthThreshold = 1.0,

        .defaultNumMemoryPorts = 2,

        .regCountMethod = PPAAnalyzerConfig::CountCarriedOverCycleBoundaries
    };

    ILA_WARN_IF(m_configuration.pipelineStartAtCycleBoundary 
        && m_configuration.pipelineShortOperations) 
        << "Setting pipelineShortOperations while pipelineStartAtCycleBoundary is set will have no effect";
    
    if (m_configuration.shortOperationLengthThreshold > 1.0)
    {
        m_configuration.shortOperationLengthThreshold = 1.0;
    }
    

    ILA_ASSERT(m_configuration.shortOperationLengthThreshold <= 1.0) 
        << "Short operation threshold is too large";

    registerAllModels(m_registrar);

    PPAAnalyze();
}

/*****************************************************************************/

PPAAnalyzer::PPAAnalyzer
(
    const InstrLvlAbsPtr & ila,
    double cycle_time,
    const std::string & instr_seq_path,
    const std::string & vcd_path,
    PPAAnalyzerConfig * config
)   : m_ila(ila), m_cycleTime(cycle_time), m_configuration(*config),
      m_instrSeqPath(instr_seq_path), m_vcdPath (vcd_path)
{


    ILA_ASSERT(m_configuration.shortOperationLengthThreshold <= 1.0) 
        << "Short operation threshold is too large";

    registerAllModels(m_registrar);

    PPAAnalyze();
}

}