/// \file
/// Implementation of DFS visitor to translate expr in Aladdin_Ilator.

#include <ilang/target-aladdin/aladdin_ilator.h>

#include <fmt/format.h>
#include <numeric>

#include <ilang/ila/ast_hub.h>
#include <ilang/util/log.h>

#include <bit>


namespace ilang {

void Aladdin_Ilator::DfsExpr(const ExprPtr& e, StrBuff& buff, ExprVarMap& lut) {
  if (auto pos = lut.find(e); pos == lut.end()) {
    if (e->is_var()) {
      DfsVar(e, buff, lut);
    } else if (e->is_const()) {
      DfsConst(e, buff, lut);
    } else {
      ILA_ASSERT(e->is_op());
      DfsOp(e, buff, lut);
    }
  }
  ILA_ASSERT((e->is_mem() && e->is_op()) || (lut.find(e) != lut.end()));
}

void Aladdin_Ilator::DfsVar(const ExprPtr& expr, StrBuff& buff,
                            ExprVarMap& lut) const {
  auto [it, status] = lut.try_emplace(expr, GetCName(expr));
  ILA_ASSERT(status);
  // no need to define new variable
}

void Aladdin_Ilator::DfsConst(const ExprPtr& expr, StrBuff& buff,
                              ExprVarMap& lut) {
  auto local_var = GetLocalVar(lut);
  auto [it, status] = lut.try_emplace(expr, local_var);
  ILA_ASSERT(status);

  // alias for constant memory
  if (expr->is_mem()) {

    const_mems_.insert(expr);

    SortPtr sort = expr->sort();

    if (sort->data_width() > 1) {
      fmt::format_to(buff,
                     "_BitInt({data_width}) * {local_var} = {const_mem};\n",
                     fmt::arg("local_var", local_var),
                     fmt::arg("data_width", sort->data_width()),
                     fmt::arg("const_mem", GetCName(expr)));
    } else {
      fmt::format_to(buff, "bool * {local_var} = {const_mem};\n",
                     fmt::arg("local_var", local_var),
                     fmt::arg("const_mem", GetCName(expr)));
    }

    // static const char* kConstMemTemplate = "{type} {local_var} =
    // {const_mem};\n"; fmt::format_to(buff, kConstMemTemplate,
    //                fmt::arg("type", GetCType(expr)),
    //                fmt::arg("local_var", local_var),
    //                fmt::arg("const_mem", GetCName(expr)));

    // auto const_mem = std::dynamic_pointer_cast<ExprConst>(expr);
    // const auto& val_map = const_mem->val_mem()->val_map();
    // static const char* kMemStoreTemplate =
    //   "_BitInt({addr_width}) {local_var} = {address};\n"
    //   "_BitInt({data_width}) {local_var}_{ind}_data = {data};\n";

    // int i = 0;
    // for (auto& it : val_map) {
    //   i++;
    //   fmt::format_to(buff, kMemStoreTemplate,
    //                   fmt::arg("local_var", local_var),
    //                   fmt::arg("address", it.first),
    //                   fmt::arg("data", it.second),
    //                   fmt::arg("addr_width",
    //                   const_mem->sort()->addr_width()),
    //                   fmt::arg("data_width",
    //                   const_mem->sort()->data_width()), fmt::arg("ind", i));
    // }

    return;
  }

  // define new var with constant value
  auto expr_const = std::dynamic_pointer_cast<ExprConst>(expr);
  std::string value = "";
  if (expr->is_bool()) {
    value = expr_const->val_bool()->val() ? "true" : "false";
  } else {
    ILA_ASSERT(expr->is_bv());
    value = std::to_string(expr_const->val_bv()->val());
  }
  static const char* kConstNonMemTemplate =
      "{var_type} {local_var} = {const_value};\n";
  fmt::format_to(buff, kConstNonMemTemplate, //
                 fmt::arg("var_type", GetCType(expr)),
                 fmt::arg("local_var", local_var),
                 fmt::arg("const_value", value));
}

void Aladdin_Ilator::DfsOp(const ExprPtr& expr, StrBuff& buff,
                           ExprVarMap& lut) {
  // store, ite-mem
  if (expr->is_mem()) {
    DfsOpMemory(expr, buff, lut);
    return;
  }

  switch (auto uid = asthub::GetUidExprOp(expr); uid) {
  // apply function
  case AstUidExprOp::kApplyFunc:
    DfsOpAppFunc(expr, buff, lut);
    break;
  // special cases
  case AstUidExprOp::kLoad:
    [[fallthrough]];
  case AstUidExprOp::kConcatenate:
    [[fallthrough]];
  case AstUidExprOp::kExtract:
    [[fallthrough]];
  case AstUidExprOp::kZeroExtend:
    [[fallthrough]];
  case AstUidExprOp::kSignedExtend:
    [[fallthrough]];
  case AstUidExprOp::kImply:
    [[fallthrough]];
  case AstUidExprOp::kIfThenElse:
    DfsOpSpecial(expr, buff, lut);
    break;
  // regular operator
  default:
    DfsOpRegular(expr, buff, lut);
    break;
  };
}

void Aladdin_Ilator::DfsOpMemory(const ExprPtr& expr, StrBuff& buff,
                                 ExprVarMap& lut) {
  // TODO CHECK DMA STUFF

  auto local_var = GetLocalVar(lut);
  auto [it, status] = lut.try_emplace(expr, local_var);
  ILA_ASSERT(status);

  if (auto uid = asthub::GetUidExprOp(expr); uid == AstUidExprOp::kStore) {

    static const char* kMemStoreTemplate =
        "{addr_type} {local_var}_addr = {address};\n"
        "{data_type} {local_var}_data = {data};\n";

    fmt::format_to(buff, kMemStoreTemplate,
                   fmt::arg("addr_type", GetCType(expr->arg(1))),
                   fmt::arg("local_var", local_var),
                   fmt::arg("address", LookUp(expr->arg(1), lut)),
                   fmt::arg("data_type", GetCType(expr->arg(2))),
                   fmt::arg("data", LookUp(expr->arg(2), lut)));

  } else { // ite

    // static const char* kMemStoreTemplate =
    //     "_BitInt({addr_width}) {local_var}_addr = {store_var}_addr;\n"
    //     "_BitInt({data_width}) {local_var}_data = {store_var}_data;\n";

    // fmt::format_to(buff, "if ({}) {{\n", LookUp(expr->arg(0), lut));
    // fmt::format_to(buff, kMemStoreTemplate,
    //                fmt::arg("addr_width", expr->sort()->addr_width()),
    //                fmt::arg("local_var", local_var),
    //                fmt::arg("store_var", LookUp(expr->arg(1), lut)),
    //                fmt::arg("data_width", expr->sort()->data_width()));
    //  fmt::format_to(buff, "}} else {{\n");
    //  fmt::format_to(buff, kMemStoreTemplate,
    //                fmt::arg("addr_width", expr->sort()->addr_width()),
    //                fmt::arg("local_var", local_var),
    //                fmt::arg("store_var", LookUp(expr->arg(2), lut)),
    //                fmt::arg("data_width", expr->sort()->data_width()));
    //   fmt::format_to(buff, "}}\n");
  }
}

void Aladdin_Ilator::DfsOpAppFunc(const ExprPtr& expr, StrBuff& buff,
                                  ExprVarMap& lut) {
  ILA_CHECK(!expr->is_mem()) << "Func returning memory not supported yet";

  auto local_var = GetLocalVar(lut);
  auto [it, status] = lut.try_emplace(expr, local_var);
  ILA_ASSERT(status);

  // apply uninterpreted function
  auto app_func = std::dynamic_pointer_cast<ExprOpAppFunc>(expr);
  auto func = app_func->func();
  auto func_cxx = RegisterExternalFunc(func);

  std::vector<std::string> arguments;
  for (size_t i = 0; i < func->arg_num(); i++) {
    arguments.push_back(LookUp(app_func->arg(i), lut));
  }

  // static const char* kAppFuncTemplate =
  //     "{return_type} {return_var} = 0;\n";
  // fmt::format_to(buff, kAppFuncTemplate, //
  //                fmt::arg("return_type", GetCType(expr)),
  //                fmt::arg("return_var", local_var));

  static const char* kAppFuncTemplate =
      "{return_type} {return_var} = {func_name}({argument_list});\n";
  fmt::format_to(buff, kAppFuncTemplate, //
                 fmt::arg("return_type", GetCType(expr)),
                 fmt::arg("return_var", local_var),
                 fmt::arg("func_name", func_cxx->name),
                 fmt::arg("argument_list", fmt::join(arguments, ", ")));
}

void Aladdin_Ilator::DfsOpSpecial(const ExprPtr& expr, StrBuff& buff,
                                  ExprVarMap& lut) {
  auto local_var = GetLocalVar(lut);
  auto [it, status] = lut.try_emplace(expr, local_var);
  ILA_ASSERT(status);

  switch (auto uid = asthub::GetUidExprOp(expr); uid) {
  case AstUidExprOp::kLoad: {
    // TODO CHECK IS DMA

    if (memory_types.find(expr->arg(0)) == memory_types.end() ||
        memory_types.at(expr->arg(0)).mt == host) {

      static const char* kLoadTemplate =
          "hostLoad(&dma_var, {memory_source} + {address}, {word_size});\n"
          "{return_type} {local_var} = ({return_type})dma_var;\n";
      uint64_t wordSize = GetWordSize(expr);
      fmt::format_to(buff, kLoadTemplate, //
                     fmt::arg("return_type", GetCType(expr)),
                     fmt::arg("local_var", local_var),
                     fmt::arg("memory_source", LookUp(expr->arg(0), lut)),
                     fmt::arg("address", LookUp(expr->arg(1), lut)),
                     fmt::arg("word_size", wordSize));
      if (wordSize > biggestDMA) {
        biggestDMA = wordSize;
      }
      if (dmaGCD == (size_t)-1) {
        dmaGCD = wordSize;
      } else {
        dmaGCD = std::gcd(dmaGCD, wordSize);
      }
    } else {
      static const char* kLoadTemplate =
          // "printf(\"%d\\n\", {address}); // to prevent it's calculation being optimized away\n"
          "{return_type} {local_var} = {memory_source}[(unsigned _BitInt({bitwidth})){address}];\n";
      fmt::format_to(buff, kLoadTemplate, //
                     fmt::arg("return_type", GetCType(expr)),
                     fmt::arg("local_var", local_var),
                     fmt::arg("memory_source", LookUp(expr->arg(0), lut)),
                    //  fmt::arg("fake_address", rand() % memory_types.at(expr->arg(0)).size),
                     fmt::arg("bitwidth", MemSizeToAddressBits(memory_types.at(expr->arg(0)).size, GetWordSize(expr->arg(0)))),
                     fmt::arg("address", LookUp(expr->arg(1), lut))

      );
    }
    break;
  }
  case AstUidExprOp::kConcatenate: {
    static const char* kConcatTemplate =
        "{return_type} {local_var} = (((({type}){arg_0}) << {arg1_width}) | (({type}){arg_1}));\n";
    auto arg0 = expr->arg(0);
    auto arg1 = expr->arg(1);
    auto arg1_width = expr->arg(1)->sort()->bit_width();
    fmt::format_to(buff, kConcatTemplate, //
                   fmt::arg("return_type", GetCType(expr)),
                   fmt::arg("local_var", local_var),
                   fmt::arg("type", GetCType(expr)),
                   fmt::arg("arg_0", LookUp(arg0, lut)),
                   fmt::arg("arg1_width", arg1_width),
                   fmt::arg("arg_1", LookUp(arg1, lut)));
    break;
  }
  // MAKE SURE THAT RETURN TYPE IS THE CORRECT WIDTH
  case AstUidExprOp::kExtract: {
    static const char* kExtractTemplate =
        "{return_type} {extract} = (({return_type}) ({origin} >> "
        "{loc_low}));\n";
    fmt::format_to(buff, kExtractTemplate, //
                   fmt::arg("return_type", GetCType(expr)),
                   fmt::arg("extract", local_var),
                   fmt::arg("origin", LookUp(expr->arg(0), lut)),
                   fmt::arg("loc_low", expr->param(1)));
    break;
  }
  // LEFT OFF HERE
  case AstUidExprOp::kZeroExtend: {
    static const char* kExtendTemplate =
        "{return_type} {extend} = (unsigned {old_type}) {origin};\n";
    auto origin_expr = expr->arg(0);
    fmt::format_to(buff, kExtendTemplate, //
                   fmt::arg("return_type", GetCType(expr)),
                   fmt::arg("extend", local_var),
                   fmt::arg("old_type", GetCType(origin_expr))),
        fmt::arg("origin", LookUp(origin_expr, lut));
    break;
  }
  case AstUidExprOp::kSignedExtend: {
    static const char* kExtendTemplate =
        "{return_type} {extend} = (signed {old_type}) {origin};\n";
    auto origin_expr = expr->arg(0);
    fmt::format_to(buff, kExtendTemplate, //
                   fmt::arg("return_type", GetCType(expr)),
                   fmt::arg("extend", local_var),
                   fmt::arg("old_type", GetCType(origin_expr))),
        fmt::arg("origin", LookUp(origin_expr, lut));
    break;
  }
  case AstUidExprOp::kImply: {
    static const char* kImplyTemplate =
        "{return_type} {local_var} = (~{if_var}) & {then_var};\n";
    fmt::format_to(buff, kImplyTemplate, //
                   fmt::arg("return_type", GetCType(expr)),
                   fmt::arg("local_var", local_var),
                   fmt::arg("if_var", LookUp(expr->arg(0), lut)),
                   fmt::arg("then_var", LookUp(expr->arg(1), lut)));
    break;
  }
  case AstUidExprOp::kIfThenElse: {
    static const char* kIteTemplate =
        "{return_type} {local_var} = ({condition}) ? {true_branch} : "
        "{false_branch};\n";
    fmt::format_to(buff, kIteTemplate, //
                   fmt::arg("return_type", GetCType(expr)),
                   fmt::arg("local_var", local_var),
                   fmt::arg("condition", LookUp(expr->arg(0), lut)),
                   fmt::arg("true_branch", LookUp(expr->arg(1), lut)),
                   fmt::arg("false_branch", LookUp(expr->arg(2), lut)));
    break;
  }
  default:
    ILA_CHECK(false) << expr;
    break;
  };
}

static const std::unordered_map<AstUidExprOp, std::string> kOpSymbols = {
    // unary
    {AstUidExprOp::kNegate, "-"},
    {AstUidExprOp::kNot, "!"},
    {AstUidExprOp::kComplement, "~"},
    // binary compare
    {AstUidExprOp::kEqual, "=="},
    {AstUidExprOp::kLessThan, "<"},
    {AstUidExprOp::kGreaterThan, ">"},
    {AstUidExprOp::kUnsignedLessThan, "<"},
    {AstUidExprOp::kUnsignedGreaterThan, ">"},
    // binary arith
    {AstUidExprOp::kAnd, "&"},
    {AstUidExprOp::kOr, "|"},
    {AstUidExprOp::kXor, "^"},
    {AstUidExprOp::kShiftLeft, "<<"},
    {AstUidExprOp::kLogicShiftRight, ">>"},
    {AstUidExprOp::kArithShiftRight, ">>"},
    {AstUidExprOp::kAdd, "+"},
    {AstUidExprOp::kSubtract, "-"},
    {AstUidExprOp::kMultiply, "*"},
    {AstUidExprOp::kDivide, "/"},
    {AstUidExprOp::kUnsignedRemainder, "%"}};

void Aladdin_Ilator::DfsOpRegular(const ExprPtr& expr, StrBuff& buff,
                                  ExprVarMap& lut) const {
  auto local_var = GetLocalVar(lut);
  auto [it, status] = lut.try_emplace(expr, local_var);
  ILA_ASSERT(status);

  // get the corresponding operator symbol
  auto uid = asthub::GetUidExprOp(expr);
  auto pos = kOpSymbols.find(uid);
  ILA_ASSERT(pos != kOpSymbols.end()) << uid;

  static const char* kUnaryOpTemplate =
      "{var_type} {local_var} = {unary_op}{arg_0};\n";
  static const char* kBinaryOpTemplate =
      "{var_type} {local_var} = ({arg_0} {binary_op} {arg_1});\n";

  if (expr->arg_num() == 1) {
    fmt::format_to(buff, kUnaryOpTemplate, //
                   fmt::arg("var_type", GetCType(expr)),
                   fmt::arg("local_var", local_var),
                   fmt::arg("unary_op", pos->second),
                   fmt::arg("arg_0", LookUp(expr->arg(0), lut)));
  } else if (expr->arg_num() == 2) {
    fmt::format_to(buff, kBinaryOpTemplate, //
                   fmt::arg("var_type", GetCType(expr)),
                   fmt::arg("local_var", local_var),
                   fmt::arg("arg_0", LookUp(expr->arg(0), lut)),
                   fmt::arg("binary_op", pos->second),
                   fmt::arg("arg_1", LookUp(expr->arg(1), lut)));
  }
  ILA_ASSERT(expr->arg_num() <= 2);
}

} // namespace ilang
