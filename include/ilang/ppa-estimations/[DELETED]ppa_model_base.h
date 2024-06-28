/// \file
///

/*
#ifndef ILANG_ILA_PPA_ESTIMATIONS_PPA_MODEL_BASE_H__
#define ILANG_ILA_PPA_ESTIMATIONS_PPA_MODEL_BASE_H__

#include <ilang/ila/instr_lvl_abs.h>
#include <ilang/ila/ast/expr_op.h>

namespace ilang {

class PPAModelBase {

public:

    virtual double getOpTime(AstUidExprOp op, double maximumTime);

    virtual double getOpDynamicPower(AstUidExprOp op);
    virtual double getOpDynamicPower(AstUidExprOp op, float switchingActivity);

    virtual double getOpLeakagePower(AstUidExprOp op);

    virtual double getOpArea(AstUidExprOp op);

};

}

#endif//ILANG_ILA_PPA_ESTIMATIONS_PPA_H__
*/