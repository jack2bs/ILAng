/// \file
/// Implementation of the class Aladdin_Ilator.

#include <ilang/target-aladdin/aladdin_ilator.h>

#include <fstream>

#include <fmt/format.h>
#include <numeric>
#include <z3++.h>

#include <ilang/config.h>
#include <ilang/ila-mngr/pass.h>
#include <ilang/ila-mngr/u_abs_knob.h>
#include <ilang/ila/ast_hub.h>
#include <ilang/target-smt/z3_expr_adapter.h>
#include <ilang/util/fs.h>
#include <ilang/util/log.h>

/// \namespace ilang
namespace ilang {

//
// static helpers/members
//

static const std::string kDirSrc = "ila";

static std::unordered_map<size_t, size_t> kPivotalId;

static size_t GetPivotalId(const size_t& id) {
  if (auto pos = kPivotalId.find(id); pos == kPivotalId.end()) {
    auto new_id = kPivotalId.size();
    kPivotalId.insert({id, new_id});
    return new_id;
  } else {
    return pos->second;
  }
}

static void WriteFile(const std::string& file_path,
                      const fmt::memory_buffer& buff) {
  std::ofstream fw(file_path);
  ILA_ASSERT(fw.is_open()) << "Fail opening file " << file_path;
  fw << to_string(buff);
  fw.close();
}

static bool HasLoadFromStore(const ExprPtr& expr) {
  auto monitor = false;
  auto LoadFromStore = [&monitor](const ExprPtr& e) {
    if (e->is_op()) {
      if (asthub::GetUidExprOp(e) == AstUidExprOp::kLoad) {
        monitor |= e->arg(0)->is_op();
      }
    }
  };
  expr->DepthFirstVisit(LoadFromStore);
  return monitor;
}

//
// Aladdin_Ilator implementation
//

Aladdin_Ilator::Aladdin_Ilator(const InstrLvlAbsPtr& m) : m_(m) {

  input_memory_type = end_mem_type;
}

Aladdin_Ilator::~Aladdin_Ilator() { Reset(); }

void Aladdin_Ilator::Generate(const std::string& dst, bool opt) {
  // sanity checks and initialize
  if (!SanityCheck() || !Bootstrap(dst, opt)) {
    return;
  }

  auto status = true;
  ILA_INFO << "Start generating Aladdin simulator of " << m_;

  // Get information from the user
  status &= GetInformationFromUser(os_portable_append_dir(dst, kDirSrc));

  // non-instruction basics
  status &= GenerateIlaBasics(os_portable_append_dir(dst, kDirSrc));

  // instruction semantics (decode and updates)
  for (auto& instr : absknob::GetInstrTree(m_)) {
    status &= GenerateInstrContent(instr, os_portable_append_dir(dst, kDirSrc));
  }

  // constant memory (Need to check)
  status &= GenerateConstantMemory(os_portable_append_dir(dst, kDirSrc));

  // initial condition setup (Need to check)
  status &= GenerateInitialSetup(os_portable_append_dir(dst, kDirSrc));

  // valid and decode simulator
  status &= GenerateSim(os_portable_append_dir(dst, kDirSrc));

  // shared header (input, state, func., etc.)
  status &= GenerateGlobalHeader(os_portable_append_dir(dst, kDirSrc));

  // file which configures the variables
  status &= GenerateConfigTemplate(os_portable_append_dir(dst, kDirSrc));

  // cmake support, e.g., recipe and templates
  status &= GenerateBuildSupport(os_portable_append_dir(dst, kDirSrc));

  // clean up if something went wrong
  if (status) {
    ILA_INFO << "Sucessfully generate Aladdin simulator at " << dst;
  } else {
    ILA_ERROR << "Fail generating simulator at " << dst;
#ifdef NDEBUG
    os_portable_remove_directory(dst);
#endif // NDEBUG
  }
}

void Aladdin_Ilator::Reset() {
  // functions
  for (auto f : functions_) {
    delete f.second;
  }
  functions_.clear();

  // externs
  for (auto f : externs_) {
    delete f.second;
  }
  externs_.clear();

  const_mems_.clear();
  sourcesList.clear();
}

bool Aladdin_Ilator::SanityCheck() const {
  // add new check here
  return true;
}

bool Aladdin_Ilator::Bootstrap(const std::string& root, bool opt) {
  Reset();
  auto status = true;

  // light-weight preprocessing
  if (opt) {
    status &= pass::SimplifySyntactic(m_);
    status &= pass::RewriteConditionalStore(m_);
  }

  // create/structure project directory
  status &= os_portable_mkdir(root);
  status &= os_portable_mkdir(os_portable_append_dir(root, kDirSrc));
  if (!status) {
    os_portable_remove_directory(root);
  }

  ILA_ERROR_IF(!status) << "Fail bootstraping";
  return status;
}

Aladdin_Ilator::MemoryType Aladdin_Ilator::GetMemoryTypeInput() {
  std::string nextLine;
  std::cin >> nextLine;

  // Make the input lowercase
  std::transform(nextLine.begin(), nextLine.end(), nextLine.begin(),
                 [](unsigned char c) { return std::tolower(c); });

  if (nextLine.find("host") != std::string::npos) {
    return host;
  }
  if (nextLine.find("spad") != std::string::npos) {
    return spad;
  }
  return end_mem_type;
}

std::string Aladdin_Ilator::MemoryTypeToString(MemoryType inp) {
  switch (inp) {
  case spad:
    return "spad";
  case host:
    return "host";
  default:
    return "end_mem_type";
  }
}

Aladdin_Ilator::MemoryType Aladdin_Ilator::StringToMemoryType(std::string inp) {

  // Make the input lowercase
  std::transform(inp.begin(), inp.end(), inp.begin(),
                 [](unsigned char c) { return std::tolower(c); });

  if (inp.find("spad") != std::string::npos) {
    return spad;
  }
  if (inp.find("host") != std::string::npos) {
    return host;
  }
  return end_mem_type;
}

bool Aladdin_Ilator::GetInformationFromUser(const std::string& dir) {

  StrBuff buff;

  char reuse = 'n';

  if (os_portable_exist(os_portable_append_dir(dir, "input_info"))) {
    std::cout << "Reuse existing choices?(y / n)\n";
    std::cin >> reuse;
  }

  if (reuse == 'y' || reuse == 'Y') {
    std::ifstream input_file(os_portable_append_dir(dir, "input_info"));
    std::string file_string = "";
    std::ostringstream ss;
    ss << input_file.rdbuf();
    file_string = ss.str();

    int pos = file_string.find("InputMemType[");
    if (pos != std::string::npos) {
      std::string asStr = file_string.substr(pos, file_string.find(']', pos));
      input_memory_type = StringToMemoryType(asStr);
      fmt::format_to(buff, "InputMemType[{}]\n",
                     MemoryTypeToString(input_memory_type));
    }

    for (int i = 0; i < m_->state_num(); i++) {
      ExprPtr var = m_->state(i);
      if (var->is_mem()) {

        int pos = file_string.find(var->name().str() + "MemType[");
        if (pos != std::string::npos) {
          std::string asStr =
              file_string.substr(pos, file_string.find(']', pos));
          memory_types.insert({var, StringToMemoryType(asStr)});
          fmt::format_to(buff, "{}MemType[{}]\n", var->name().str(),
                         MemoryTypeToString(memory_types.at(var)));
        }
      }
    }
  }

  if (input_memory_type == end_mem_type) {
    std::cout << "What memory type should be used for the inputs?"
                 " (host / spad)\n";

    input_memory_type = GetMemoryTypeInput();
    if (input_memory_type == end_mem_type) {
      std::cerr << "Invalid memory type given for inputs\n";
      exit(-1);
    }
    fmt::format_to(buff, "InputMemType[{}]\n",
                   MemoryTypeToString(input_memory_type));
  }
  for (int i = 0; i < m_->state_num(); i++) {
    ExprPtr var = m_->state(i);
    {
      if (var->is_mem()) {
        if (memory_types.find(var) == memory_types.end() ||
            memory_types.at(var) == end_mem_type) {
          std::cout << "What memory type should be used for state \""
                    << var->name().str() << "\" (host / spad)\n";
          memory_types.insert({var, GetMemoryTypeInput()});

          if (memory_types.find(var) != memory_types.end() &&
              memory_types.at(var) == end_mem_type) {
            std::cerr << "Invalid memory type given for state \""
                      << var->name().str() << "\"\n";
            exit(-1);
          }
          fmt::format_to(buff, "{}MemType[{}]\n", var->name().str(),
                         MemoryTypeToString(memory_types.at(var)));
        }
      }
    }
  }

  CommitSource("input_info", dir, buff);

  return true;
}

bool Aladdin_Ilator::GenerateIlaBasics(const std::string& dir) {
  StrBuff buff;

  // include headers
  fmt::format_to(buff, "#include \"{}.h\"\n", GetProjectName());

  // generate valid func for each ILA
  auto PerIla = [this, &buff](const InstrLvlAbsCnstPtr& m) {
    ILA_NOT_NULL(m);
    auto valid_expr = m->valid();
    if (!valid_expr) {
      valid_expr = asthub::BoolConst(true);
      ILA_WARN << "Use default (true) valid for " << m;
    }
    static int valid_count = 0;
    ExprVarMap lut;

    // Render it once
    StrBuff rendrBuff;
    ILA_CHECK(RenderExpr(valid_expr, rendrBuff, lut));

    // Put it in the valid and decode paralellized function
    fmt::format_to(validAndDecodeFunc, "bool valid{} = false;\n", valid_count);
    fmt::format_to(validAndDecodeFunc, "{{\n");
    fmt::format_to(validAndDecodeFunc, "{}", fmt::to_string(rendrBuff));
    fmt::format_to(validAndDecodeFunc, "valid{} = {};\n", valid_count,
                   LookUp(valid_expr, lut));
    fmt::format_to(validAndDecodeFunc, "}}\n");
    valid_count++;

    // Put it in it's own function
    auto valid_func = RegisterFunction(GetValidFuncName(m), valid_expr);
    BeginFuncDef(valid_func, buff);

  fmt:;
    format_to(buff, "{}", fmt::to_string(rendrBuff));
    fmt::format_to(buff, "{var} = {local_name};\n",
                   fmt::arg("var", GetCNameWithType(valid_expr)),
                   fmt::arg("local_name", LookUp(valid_expr, lut)));
    EndFuncDef(valid_func, buff);
  };

  // traverse the hierarchy
  m_->DepthFirstVisit(PerIla);

  // record and write the individual valid functions to a file
  CommitSource("all_valid_funcs_in_hier.c", dir, buff);
  return true;
}

void Aladdin_Ilator::DMAVarAdj(int wordSize) {
  if (wordSize > biggestDMA) {
    biggestDMA = wordSize;
  }

  if (dmaGCD == (size_t)-1) {
    dmaGCD = wordSize;
  } else {
    dmaGCD = std::gcd(dmaGCD, wordSize);
  }
}

bool Aladdin_Ilator::MemHelper(StrBuff& b, ExprVarMap& l, ExprPtr& old,
                               ExprPtr& next) {

  // helper for traversing memory updates
  class MemUpdateVisiter {
  public:
    MemUpdateVisiter(Aladdin_Ilator* h, StrBuff& b, ExprVarMap& l, ExprPtr& r,
                     bool d)
        : host(h), buff(b), lut_ref(l), root(r), isHostMem(d) {}

    bool pre(const ExprPtr& expr) {
      // stop traversing when reaching memory ITE (stand-alone func)
      if (expr->is_mem() && expr->is_op() &&
          asthub::GetUidExprOp(expr) == AstUidExprOp::kIfThenElse) {

        auto pos = lut_ref.find(expr->arg(0));
        ILA_ASSERT(pos != lut_ref.end());
        std::string cond = pos->second;

        fmt::format_to(buff, "if ({}){{\\\n", cond);
        expr->arg(1)->DepthFirstVisitPrePost(*this);
        fmt::format_to(buff, "}} else {{\\\n");
        expr->arg(2)->DepthFirstVisitPrePost(*this);
        fmt::format_to(buff, "}}\\\n");
        // host->DfsExpr(expr, buff_ref, lut_ref);
        return true;
      } else {
        return false;
      }
    }

    void post(const ExprPtr& expr) {
      if (expr->is_mem() && expr->is_op() &&
          asthub::GetUidExprOp(expr) == AstUidExprOp::kStore) {

        if (isHostMem) {
          auto pos = lut_ref.find(expr);
          ILA_ASSERT(pos != lut_ref.end());
          uint64_t wordSize = GetWordSize(root);
          std::string local_var = pos->second;
          fmt::format_to(buff,
                         "dma_var = {local_var}_data;\n"
                         "hostStore({mem_name} + {local_var}_addr, &dma_var, "
                         "{word_size});\n",
                         fmt::arg("mem_name", GetCName(root)),
                         fmt::arg("local_var", local_var),
                         fmt::arg("word_size", wordSize));

          host->DMAVarAdj(wordSize);

        } else {
          auto pos = lut_ref.find(expr);
          ILA_ASSERT(pos != lut_ref.end());
          std::string local_var = pos->second;
          fmt::format_to(buff,
                         "{mem_name}[{local_var}_addr] = {local_var}_data;\n",
                         fmt::arg("mem_name", GetCName(root)),
                         fmt::arg("local_var", local_var));
        }

      } else if (expr->is_mem() && expr->is_const()) {
        auto const_mem = std::dynamic_pointer_cast<ExprConst>(expr);
        const auto& val_map = const_mem->val_mem()->val_map();

        if (isHostMem) {
          // auto pos = lut_ref.find(expr);
          // ILA_ASSERT(pos != lut_ref.end());
          uint64_t wordSize = GetWordSize(root);
          // std::string local_var = pos->second;
          for (auto& it : val_map) {
            fmt::format_to(
                buff,
                "dma_var = {data};\n"
                "hostStore({mem_name} + {addr}, &dma_var, {word_size});\n",
                fmt::arg("mem_name", GetCName(root)),
                fmt::arg("word_size", wordSize), fmt::arg("data", it.second),
                fmt::arg("addr", it.first));
          }

          host->DMAVarAdj(wordSize);

        } else {
          // auto pos = lut_ref.find(expr);
          // ILA_ASSERT(pos != lut_ref.end());
          // std::string local_var = pos->second;
          for (auto& it : val_map) {

            fmt::format_to(buff, "{mem_name}[{addr}] = {data};\n",
                           fmt::arg("mem_name", GetCName(root)),
                           fmt::arg("addr", it.first),
                           fmt::arg("data", it.second));
          }
        }
      }
    }

  private:
    Aladdin_Ilator* host;
    StrBuff& buff;
    ExprVarMap lut_ref;
    ExprPtr& root;
    bool isHostMem;
  };

  auto RenderMemUpdate = [this](const ExprPtr& e, StrBuff& b, ExprVarMap& l,
                                ExprPtr& r) {
    bool d = memory_types.find(r) != memory_types.end() &&
             memory_types.at(r) == host;

    auto mem_visiter = MemUpdateVisiter(this, b, l, r, d);
    e->DepthFirstVisitPrePost(mem_visiter);
  };

  RenderMemUpdate(next, b, l, old);

  return true;
}

bool Aladdin_Ilator::GenerateInstrContent(const InstrPtr& instr,
                                          const std::string& dir) {
  StrBuff buff;

  ExprVarMap lut;
  StrBuff decodeRendrBuff;
  auto decode_expr = instr->decode();
  if (!RenderExpr(decode_expr, decodeRendrBuff, lut)) {
    return false;
  }

  // Put the decode expression in the parallelized valid and decode function
  static int instr_count = 0;
  fmt::format_to(validAndDecodeFunc, "bool decode{} = false;\n", instr_count);
  fmt::format_to(validAndDecodeFunc, "{{\n");
  fmt::format_to(validAndDecodeFunc, "{}", fmt::to_string(decodeRendrBuff));
  fmt::format_to(validAndDecodeFunc, "decode{} = {};\n", instr_count,
                 LookUp(decode_expr, lut));
  fmt::format_to(validAndDecodeFunc, "}}\n");
  instr_count++;

  fmt::format_to(buff, "#include \"{}.h\"\n", GetProjectName());

  // Put the decode in it's own function
  auto decode_func = RegisterFunction(GetDecodeFuncName(instr), decode_expr);
  BeginFuncDef(decode_func, buff);
  fmt::format_to(buff, "{}", fmt::to_string(decodeRendrBuff));
  fmt::format_to(buff, "{var} = {local_name};\n",
                 fmt::arg("var", GetCNameWithType(decode_expr)),
                 fmt::arg("local_name", LookUp(decode_expr, lut)));
  EndFuncDef(decode_func, buff);

  //
  // next state

  lut.clear();
  StrBuff updateFuncBuff;
  std::set<ExprPtr> visited;
  auto updated_states = instr->updated_states();
  for (const auto& s : updated_states) {
    // check if visited
    auto update_expr = instr->update(s);
    if (auto pos = visited.find(update_expr); pos == visited.end()) {
      visited.insert(update_expr);
    } else {
      continue;
    }
    // create placeholder
    if (auto update_expr = instr->update(s); !update_expr->is_mem()) {
      if (!RenderExpr(update_expr, updateFuncBuff, lut)) {
        return false;
      }
      fmt::format_to(updateFuncBuff,
                     "{type_name} {local_var}_nxt_holder = {local_var};\n",
                     fmt::arg("type_name", GetCType(update_expr)),
                     fmt::arg("local_var", LookUp(update_expr, lut)));
    } else {
      if (HasLoadFromStore(update_expr)) {
        // Load from store is not allowed
        return false;
      }
      if (!RenderExpr(update_expr, updateFuncBuff, lut)) {
        return false;
      }
    }
  }
  // update state
  for (auto& s : updated_states) {
    auto curr = instr->host()->state(s);
    auto next = instr->update(s);
    if (!curr->is_mem()) {
      fmt::format_to(updateFuncBuff, "{current} = {next_value}_nxt_holder;\n",
                     fmt::arg("current", GetCName(curr)),
                     fmt::arg("next_value", LookUp(next, lut)));
    } else {

      MemHelper(updateFuncBuff, lut, curr, next);
      continue;
    }
  }

  auto update_func = RegisterFunction(GetUpdateFuncName(instr));
  BeginFuncDef(update_func, buff);
  fmt::format_to(buff, "{}", fmt::to_string(updateFuncBuff));
  EndFuncDef(update_func, buff);

  auto ex = new ilang::ExprConst(true);
  auto ptrEx = ilang::ExprPtr(ex);
  auto update_func_seq =
      RegisterFunction(GetUpdateFuncName(instr) + "_seq", ptrEx);
  BeginFuncDef(update_func_seq, buff);
  fmt::format_to(buff, "{}", fmt::to_string(updateFuncBuff));
  fmt::format_to(buff, "return true;\n}}\n");
  // EndFuncDef(update_func_seq, buff);

  // record and write to file
  CommitSource(fmt::format("idu_{}.c", instr->name().str()), dir, buff);
  return true;
}

bool Aladdin_Ilator::GenerateConstantMemory(const std::string& dir) {
  StrBuff buff;
  fmt::format_to(buff, "#include \"{}.h\"\n", GetProjectName());

  for (auto& mem : const_mems_) {
    auto const_mem = std::dynamic_pointer_cast<ExprConst>(mem);
    const auto& val_map = const_mem->val_mem()->val_map();
    std::vector<std::string> addr_data_pairs;

    for (auto& it : val_map) {
      addr_data_pairs.push_back(fmt::format(
          "{mem_name}[{addr}] = {data};", fmt::arg("mem_name", GetCName(mem)),
          fmt::arg("addr", it.first), fmt::arg("data", it.second)));
    }
    fmt::format_to(
        buff,
        "{var};\n"
        "{addr_data_pairs}\n",
        fmt::arg("var", GetCNameWithType(mem)),
        fmt::arg("addr_data_pairs", fmt::join(addr_data_pairs, "\n")));
  }

  CommitSource("constant_memory_def.c", dir, buff);
  return true;
}

bool Aladdin_Ilator::GenerateInitialSetup(const std::string& dir) {
  // conjunct all initial condition
  auto init = asthub::BoolConst(true);
  auto ConjInit = [&init](const InstrLvlAbsCnstPtr& m) {
    for (size_t i = 0; i < m->init_num(); i++) {
      init = asthub::And(init, m->init(i));
    }
  };
  m_->DepthFirstVisit(ConjInit);

  // get value for referred vars
  z3::context ctx;
  z3::solver solver(ctx);
  Z3ExprAdapter gen(ctx);
  solver.add(gen.GetExpr(init));
  auto res = solver.check();
  if (res != z3::sat) {
    ILA_ERROR << "Fail finding assignment satisfying initial condition";
    return false;
  }

  std::map<ExprPtr, uint64_t> init_values;
  auto model = solver.get_model();
  auto refer_vars = absknob::GetVar(init);
  for (const auto& var : refer_vars) {
    auto var_value = model.eval(gen.GetExpr(var));
    try {
#ifndef Z3_LEGACY_API
      auto value_holder = var_value.get_numeral_uint64();
#else
      __uint64 value_holder;
      Z3_get_numeral_uint64(ctx, var_value, &value_holder);
#endif
      init_values.emplace(var, value_holder);
    } catch (...) {
      ILA_ERROR << "Fail getting " << var_value;
      return false;
    }
  }

  // gen file
  auto init_func = RegisterMacro("setup_initial_condition");
  StrBuff buff;

  BeginMacroDef(init_func, buff);
  for (auto pair : init_values) {
    fmt::format_to(buff, "{var_name} = {var_value};\\\n",
                   fmt::arg("var_name", GetCName(pair.first)),
                   fmt::arg("var_value", pair.second));
  }
  EndMacroDef(init_func, buff);
  CommitSource("setup_initial_condition.h", dir, buff);
  return true;
}

void Aladdin_Ilator::AddParametersToBuffer(StrBuff& b, bool withDataTypes,
                                           std::string separator) const {

  bool first = true;
  // states
  for (auto& var : absknob::GetSttTree(m_)) {
    if (var->is_mem()) {
      if (!first) {
        fmt::format_to(b, separator);
      } else {
        first = false;
      }
      if (var->sort()->data_width() > 1) {
        if (withDataTypes) {
          fmt::format_to(b, "_BitInt({data_width}) * ",
                         fmt::arg("data_width", var->sort()->data_width()));
        }

        fmt::format_to(b, "{name}", fmt::arg("name", GetCName(var)));
      } else {
        if (withDataTypes) {
          fmt::format_to(b, "bool * ");
        }

        fmt::format_to(b, "{name}", fmt::arg("name", GetCName(var)));
      }
    }
  }

  // inputs
  for (size_t i = 0; i < m_->input_num(); i++) {
    if (!first) {
      fmt::format_to(b, separator);
    } else {
      first = false;
    }

    fmt::format_to(b, "{var}",
                   fmt::arg("var", withDataTypes
                                       ? GetCNameWithType(m_->input(i))
                                       : GetCName(m_->input(i))));
  }

  fmt::format_to(b, "\n");
}

bool Aladdin_Ilator::GenerateSim(const std::string& dir) {
  StrBuff buff;

  fmt::format_to(buff, "#include \"{}.h\"\n", GetProjectName());
  fmt::format_to(buff, "#include \"gem5/dma_interface.h\"\n");
  fmt::format_to(buff, "#include <stdio.h>\n");
  fmt::format_to(buff, "#include <stdlib.h>\n");
  fmt::format_to(buff, "#define ARGSLIST \\\n");

  AddParametersToBuffer(buff, false, ",\\\n");

  // states
  if (absknob::GetSttTree(m_).size()) {
    fmt::format_to(buff, "\n// states\n");
  }
  for (auto& var : absknob::GetSttTree(m_)) {
    if (!var->is_mem()) {
      fmt::format_to(buff, "{var} = 1;\n",
                     fmt::arg("var", GetCNameWithType(var)));
    }
  }

  fmt::format_to(buff, "\n_BitInt({biggestDma}) dma_var;\n",
                 fmt::arg("biggestDma", biggestDMA * 8));

  fmt::format_to(buff, "\nbool valid_decode_sim(");

  AddParametersToBuffer(buff, true, ",\n  ");

  fmt::format_to(buff, ")");

  fmt::format_to(buff, "{{\n");

  int valNum = 0;

  // Valid for each ILA
  auto PerIla = [this, &buff, &valNum](const InstrLvlAbsCnstPtr& m) {
    fmt::format_to(buff,
                   "  bool valid{num} = false;\n"
                   "  valid{num} = {valid_func_name}(ARGSLIST);\n\n",
                   fmt::arg("num", valNum),
                   fmt::arg("valid_func_name", GetValidFuncName(m)));
    valNum++;
  };

  auto top_instrs = absknob::GetInstr(m_);
  auto all_instrs = absknob::GetInstrTree(m_);

  m_->DepthFirstVisit(PerIla);

  int decNum = 0;

  // Decode for each instruction
  auto DecInstr = [this, &buff](const InstrPtr& instr, bool child, int num) {
    fmt::format_to(buff,
                   "  bool decode{num} = false;\n"
                   "  decode{num} = {decode_func_name}(ARGSLIST);\n\n",
                   fmt::arg("num", num),
                   fmt::arg("decode_func_name", GetDecodeFuncName(instr)));
  };

  for (auto& instr : all_instrs) {
    DecInstr(instr, false, decNum);
    decNum++;
  }

  bool first = true;
  fmt::format_to(buff, "  return(\n");
  for (int i = 0; i < valNum; i++) {
    if (!first) {
      fmt::format_to(buff, "^\n");
    } else {
      first = false;
    }
    fmt::format_to(buff, "  valid{num}", fmt::arg("num", i));
  }
  for (int i = 0; i < decNum; i++) {
    if (!first) {
      fmt::format_to(buff, "^\n");
    } else {
      first = false;
    }
    fmt::format_to(buff, "  decode{num}", fmt::arg("num", i));
  }

  fmt::format_to(buff, "\n  );\n}}\n");

  fmt::format_to(buff, "\nbool valid_decode_sim_parallel(");
  AddParametersToBuffer(buff, true, ",\n  ");
  fmt::format_to(buff, ")");
  fmt::format_to(buff, "{{\n");
  fmt::format_to(buff, "{}", fmt::to_string(validAndDecodeFunc));
  first = true;
  fmt::format_to(buff, "  return(\n");
  for (int i = 0; i < valNum; i++) {
    if (!first) {
      fmt::format_to(buff, "^\n");
    } else {
      first = false;
    }
    fmt::format_to(buff, "  valid{num}", fmt::arg("num", i));
  }
  for (int i = 0; i < decNum; i++) {
    if (!first) {
      fmt::format_to(buff, "^\n");
    } else {
      first = false;
    }
    fmt::format_to(buff, "  decode{num}", fmt::arg("num", i));
  }

  fmt::format_to(buff, "\n  );\n}}\n");

  fmt::format_to(buff, "\nvoid get_base_leakage()\n"
                       "{{\n"
                       "  printf(\"fn for base leakage\\n\");\n"
                       "}}\n");

  fmt::format_to(buff, "\nvoid update_sim(");
  AddParametersToBuffer(buff, true, ",\n  ");
  fmt::format_to(buff, ")");
  fmt::format_to(buff, "{{\n");
  auto UpdInstr = [this, &buff](const InstrPtr& instr, bool child, int num) {
    fmt::format_to(buff, "  {update_func_name}(ARGSLIST);\n\n",
                   fmt::arg("update_func_name", GetUpdateFuncName(instr)));
  };
  for (auto& instr : all_instrs) {
    UpdInstr(instr, false, decNum);
  }
  fmt::format_to(buff, "\n}}\n");

  fmt::format_to(buff, "\nbool update_sim_seq(");
  AddParametersToBuffer(buff, true, ",\n  ");
  fmt::format_to(buff, ")");
  fmt::format_to(buff, "{{\n");
  auto UpdInstrSeq = [this, &buff](const InstrPtr& instr, bool child, int num) {
    fmt::format_to(
        buff, "  {seq_update_func_name}(ARGSLIST);\n\n",
        fmt::arg("seq_update_func_name", GetUpdateFuncName(instr) + "_seq"));
  };
  for (auto& instr : all_instrs) {
    UpdInstrSeq(instr, false, decNum);
  }
  fmt::format_to(buff, "\n}}\n");

  fmt::format_to(buff, "\nvoid sim(");
  AddParametersToBuffer(buff, true, ",\n  ");
  fmt::format_to(buff, ")");
  fmt::format_to(buff, "{{\n");
  fmt::format_to(buff, "  volatile bool retVal = valid_decode_sim(ARGSLIST);\n"
                       "  retVal = valid_decode_sim_parallel(ARGSLIST);\n"
                       "  update_sim(ARGSLIST);\n"
                       "  setup_initial_condition();\n"
                       "  update_sim_seq(ARGSLIST);\n"
                       "  get_base_leakage();\n");
  fmt::format_to(buff, "}}\n");

  StrBuff memoryAllocateBuff;
  memoryAllocateBuff.clear();

  // state
  if (absknob::GetSttTree(m_).size()) {
    fmt::format_to(memoryAllocateBuff, "\n  // states\n");
  }
  for (auto& var : absknob::GetSttTree(m_)) {
    // fmt::format_to(memoryAllocateBuff, "{var};\n",
    //                fmt::arg("var", GetCNameWithType(var)));
    if (var->is_mem()) {
      if (var->sort()->data_width() > 1) {
        fmt::format_to(memoryAllocateBuff,
                       "  _BitInt({width}) * {var} = (_BitInt({width})*) "
                       "malloc({size});\n",
                       fmt::arg("var", GetCName(var)),
                       fmt::arg("width", var->sort()->data_width()),
                       fmt::arg("size", GetNumBytes(var)));
      } else {
        fmt::format_to(
            memoryAllocateBuff, "  bool * {var} = (bool*) malloc({size});\n",
            fmt::arg("var", GetCName(var)), fmt::arg("size", GetNumBytes(var)));
      }
    }
  }

  // memory constant
  if (const_mems_.size()) {
    fmt::format_to(memoryAllocateBuff, "\n  // memory constants\n");
  }

  for (auto& mem : const_mems_) {
    auto const_mem = std::dynamic_pointer_cast<ExprConst>(mem);
    const auto& val_map = const_mem->val_mem()->val_map();
    std::vector<std::string> addr_data_pairs;

    for (auto& it : val_map) {
      addr_data_pairs.push_back(fmt::format(
          "{mem_name}[{addr}] = {data};", fmt::arg("mem_name", GetCName(mem)),
          fmt::arg("addr", it.first), fmt::arg("data", it.second)));
    }
    fmt::format_to(
        memoryAllocateBuff,
        "  {var};\n"
        "  {addr_data_pairs}\n",
        fmt::arg("var", GetCNameWithType(mem)),
        fmt::arg("addr_data_pairs", fmt::join(addr_data_pairs, "\n  ")));
  }

  StrBuff argsBuff;
  first = true;
  // states
  for (auto& var : absknob::GetSttTree(m_)) {
    if (var->is_mem()) {
      if (!first) {
        fmt::format_to(argsBuff, ",\n");
      } else {
        first = false;
      }

      fmt::format_to(argsBuff, "  {name}", fmt::arg("name", GetCName(var)));
    }
  }
  for (size_t i = 0; i < m_->input_num(); i++) {
    if (!first) {
      fmt::format_to(argsBuff, ",\n");
    } else {
      first = false;
    }
    fmt::format_to(argsBuff, "  rand()",
                   fmt::arg("var", GetCNameWithType(m_->input(i))));
  }

  // main function
  static const char* kSimEntryTemplate =
      "\nint main(int argc, char* argv[]) {{\n"
      "{memory_allocate}\n"
      "  setup_initial_condition();\n"
      "#ifdef LLVM_TRACE\n"
      "  sim({compute_args});\n"
      "#else\n"
      "{array_mapping}"
      "  invokeAcceleratorAndBlock(INTEGRATION_TEST);\n"
      "#endif\n"
      "  return 0; \n"
      "}}\n";

  StrBuff mappingArrays;
  for (auto s : memory_types) {
    if (s.second != MemoryType::spad) {
      uint64_t num_bytes = GetNumBytes(s.first);
      fmt::format_to(
          mappingArrays,
          "  mapArrayToAccelerator(\n"
          "    INTEGRATION_TEST, \"{var_name}\", {var_name}, {num_bytes});\n",
          fmt::arg("var_name", GetCName(s.first)),
          fmt::arg("num_bytes", num_bytes));
    }
  }

  fmt::format_to(buff, kSimEntryTemplate,
                 fmt::arg("array_mapping", to_string(mappingArrays)),
                 fmt::arg("memory_allocate", to_string(memoryAllocateBuff)),
                 fmt::arg("compute_args", to_string(argsBuff)));

  // done
  // EndFuncDef(kernel_func, buff);

  CommitSource("sim.c", dir, buff);
  return true;
}

bool Aladdin_Ilator::GenerateGlobalHeader(const std::string& dir) {
  StrBuff buff;

  fmt::format_to(buff,
                 "#ifndef _{project}_H_\n"
                 "#define _{project}_H_\n",
                 fmt::arg("project", GetProjectName()));

  fmt::format_to(buff, "\n#include <stdio.h>\n");
  fmt::format_to(buff, "\n#include <stdbool.h>\n");
  fmt::format_to(buff, "#include \"gem5/dma_interface.h\"\n");
  // fmt::format_to(buff, "#include \"gem5_harness.h\"\n");
  fmt::format_to(buff, "#include \"setup_initial_condition.h\"\n");

  fmt::format_to(buff, "\nextern _BitInt({biggestDma}) dma_var;\n",
                 fmt::arg("biggestDma", biggestDMA * 8));

  // Global state variables
  fmt::format_to(buff, "\n// states\n");
  for (auto& var : absknob::GetSttTree(m_)) {

    if (!var->is_mem()) {
      fmt::format_to(buff, "extern {var};\n",
                     fmt::arg("var", GetCNameWithType(var)));
    }
  }

  // memory constants
  fmt::format_to(buff, "\n// memory constants\n");
  for (auto& mem : const_mems_) {
    fmt::format_to(buff, "extern {var};\n",
                   fmt::arg("var", GetCNameWithType(mem)));
  }

  // function declarations
  fmt::format_to(buff, "\n// function declarations\n");
  for (auto& func : functions_) {
    WriteFuncDecl(func.second, buff);
  }
  for (auto& func : externs_) {
    WriteFuncDecl(func.second, buff);
  }

  fmt::format_to(buff, "\n#endif//_{project}_H_\n",
                 fmt::arg("project", GetProjectName()));

  // write to file
  auto file_path = os_portable_append_dir(dir, GetProjectName() + ".h");
  WriteFile(file_path, buff);
  return true;
}

bool Aladdin_Ilator::GenerateBuildSupport(const std::string& dir) {

  StrBuff buff;
  buff.clear();

  StrBuff workload;
  std::string workload_template = "{fn_name},$\\\n";

  auto all_instrs = absknob::GetInstrTree(m_);
  for (auto& instr : all_instrs) {
    fmt::format_to(workload, workload_template,
                   fmt::arg("fn_name", GetUpdateFuncName(instr)));
    fmt::format_to(workload, workload_template,
                   fmt::arg("fn_name", GetDecodeFuncName(instr)));
  }

  auto PerIla = [this, &workload](const InstrLvlAbsCnstPtr& m) {
    std::string workload_template = "{fn_name},$\\\n";

    fmt::format_to(workload, workload_template,
                   fmt::arg("fn_name", GetValidFuncName(m)));
  };
  fmt::format_to(workload, workload_template,
                 fmt::arg("fn_name", "get_base_leakage"));

  fmt::format_to(workload, workload_template,
                 fmt::arg("fn_name", "valid_decode_sim_parallel"));

  fmt::format_to(workload, workload_template,
                 fmt::arg("fn_name", "update_sim_seq"));

  m_->DepthFirstVisit(PerIla);

  fmt::format_to(buff,
                 "SRCS={src_text}\n\n"
                 "ACCEL_NAME = ila\n"
                 "TEST_BIN = $(ACCEL_NAME)\n"
                 "BMARK_SPECIFIC_CFLAGS=-DDMA_INTERFACE_V3\n"
                 "export TRACE_OUTPUT_DIR=$(ACCEL_NAME)\n"
                 "ifndef WORKLOAD\n"
                 "  export WORKLOAD={workload}\n"
                 "endif\n"
                 "include Makefile.tracer\n",
                 fmt::arg("src_text", sourcesList),
                 fmt::arg("workload", fmt::to_string(workload)));

  WriteFile(os_portable_append_dir(dir, "Makefile"), buff);

  std::string file_path = __FILE__;
  std::string dir_path =
      file_path.substr(0, file_path.rfind("aladdin_ilator.cc") - 1);
  os_portable_copy_file_to_dir(
      os_portable_append_dir(dir_path, "Makefile.tracer"), dir);

  return true;
}

void Aladdin_Ilator::AddConfigLineToBuff(const ilang::ExprPtr& var,
                                         StrBuff& buff) {
  if (!var->is_mem()) {
    fmt::format_to(buff, "partition,complete,{var},{num_bytes},{word_size}\n",
                   fmt::arg("var", GetCName(var)),
                   fmt::arg("num_bytes", GetNumBytes(var)),
                   fmt::arg("word_size", GetWordSize(var)));
  } else {
    if (memory_types.find(var) == memory_types.end()) {
      fmt::format_to(buff, "partition,cyclic,{var},{word_size},{word_size},1\n",
                     fmt::arg("var", GetCName(var)),
                     fmt::arg("word_size", GetWordSize(var)));
      return;
    }

    switch (memory_types.at(var)) {
    case spad:
      fmt::format_to(buff, "partition,cyclic,{var},{num_bytes},{word_size},1\n",
                     fmt::arg("var", GetCName(var)),
                     fmt::arg("num_bytes", GetNumBytes(var)),
                     fmt::arg("word_size", GetWordSize(var)));
      break;
    }
  }
}

bool Aladdin_Ilator::GenerateConfigTemplate(const std::string& dir) {

  StrBuff buff;

  // inputs
  // fmt::format_to(buff, "\n// inputs\n");
  // for (auto& var : absknob::GetInp(m_)) {
  //   switch (input_memory_type)
  //   {
  //   case cache:
  //     // AddConfigLineToBuff(var, buff);
  //     fmt::format_to(buff,
  //                "cache,{var}_inps_vals,{word_size}\n",
  //                fmt::arg("var", GetCName(var)),
  //                fmt::arg("num_bytes", GetNumBytes(var)),
  //                fmt::arg("word_size", GetWordSize(var)));
  //     break;
  //   case acp:
  //     // AddConfigLineToBuff(var, buff);
  //     fmt::format_to(buff,
  //                "acp,{var}_inps_vals,{word_size}\n",
  //                fmt::arg("var", GetCName(var)),
  //                fmt::arg("num_bytes", GetNumBytes(var)),
  //                fmt::arg("word_size", GetWordSize(var)));
  //     break;
  //   case spad:
  //     fmt::format_to(buff,
  //                "partition,cyclic,{var}_inps_vals,{num_bytes},{word_size},1\n",
  //                fmt::arg("var", GetCName(var)),
  //                fmt::arg("num_bytes", GetWordSize(var) * 5),
  //                fmt::arg("word_size", GetWordSize(var)));
  //     break;
  //   case dma:
  //     fmt::format_to(buff,
  //                "partition,complete,{var},{word_size},{word_size}\n",
  //                fmt::arg("var", GetCName(var)),
  //                fmt::arg("num_bytes", GetNumBytes(var)),
  //                fmt::arg("word_size", GetWordSize(var)));
  //     break;
  //   default:
  //     break;
  //   }

  // }

  // state and global vars
  // fmt::format_to(buff, "\n// states\n");
  for (auto& var : absknob::GetSttTree(m_)) {
    AddConfigLineToBuff(var, buff);
  }
  // fmt::format_to(buff, "\n// global vars\n");
  // for (auto& var : global_vars_) {
  //   AddConfigLineToBuff(var, buff);
  // }
  // fmt::format_to(buff,
  //                "partition,complete,cycles_to_simulate,4,4\n");

  // memory constant
  // fmt::format_to(buff, "\n// memory constants\n");
  for (auto& mem : const_mems_) {
    AddConfigLineToBuff(mem, buff);
  }

  fmt::format_to(buff, "partition,complete,dma_var,{},{}\n", biggestDMA,
                 dmaGCD);

  // fmt::format_to(buff,
  //                 "cycle_time,6\n"
  //                 "unrolling,compute,main_while,1\n");

  fmt::format_to(buff, "cycle_time,6\n");

  auto entry_path = os_portable_append_dir(dir, "config");

  WriteFile(entry_path, buff);
  return true;
}

bool Aladdin_Ilator::RenderExpr(const ExprPtr& expr, StrBuff& buff,
                                ExprVarMap& lut) {

  class ExprDfsVisiter {
  public:
    ExprDfsVisiter(Aladdin_Ilator* hi, StrBuff& bi, ExprVarMap& li)
        : host(hi), b(bi), l(li) {}
    bool pre(const ExprPtr& e) { return l.find(e) != l.end(); }
    void post(const ExprPtr& e) { host->DfsExpr(e, b, l); }

    Aladdin_Ilator* host;
    StrBuff& b;
    ExprVarMap& l;
  };

  try {
    auto visiter = ExprDfsVisiter(this, buff, lut);
    expr->DepthFirstVisitPrePost(visiter);
  } catch (std::exception& err) {
    ILA_ERROR << err.what();
    return false;
  }
  return true;
}

Aladdin_Ilator::CFunc*
Aladdin_Ilator::RegisterFunction(const std::string& func_name,
                                 ExprPtr return_expr) {
  auto func = new CFunc(func_name, return_expr);
  auto [it, status] = functions_.insert({func->name, func});
  ILA_ASSERT(status);
  return func;
}

Aladdin_Ilator::CFunc*
Aladdin_Ilator::RegisterMacro(const std::string& macro_name,
                              ExprPtr return_expr) {
  auto macro = new CFunc(macro_name, return_expr);
  auto [it, status] = macros_.insert({macro->name, macro});
  ILA_ASSERT(status);
  return macro;
}

Aladdin_Ilator::CFunc*
Aladdin_Ilator::RegisterExternalFunc(const FuncPtr& func) {
  auto func_cxx = new CFunc(func->name().str(), func->out());
  auto [it, status] = externs_.insert({func_cxx->name, func_cxx});

  static bool anyExterns = false;
  if (!anyExterns) {
    sourcesList += " uninterpreted_func.c";
    anyExterns = true;
  }

  // uninterpreted function can have multiple occurrence
  if (status) {
    for (size_t i = 0; i < func->arg_num(); i++) {
      it->second->args.push_back(func->arg(i));
    }
  } else {
    delete func_cxx;
  }
  return it->second;
}

void Aladdin_Ilator::BeginFuncDef(Aladdin_Ilator::CFunc* func,
                                  StrBuff& buff) const {
  ILA_ASSERT(func->args.empty()); // no definition for uninterpreted funcs

  StrBuff paramBuff;
  AddParametersToBuffer(paramBuff, true, ",\n  ");

  auto type = (func->ret) ? GetCType(func->ret) : GetCType(func->ret_type);
  auto args = (func->target)
                  ? fmt::format("{} tmp_memory", GetCType(func->target))
                  : "";

  fmt::format_to(buff, "{return_type} {func_name}({argument}) {{\n",
                 fmt::arg("return_type", type),
                 fmt::arg("project", GetProjectName()),
                 fmt::arg("func_name", func->name),
                 fmt::arg("argument", fmt::to_string(paramBuff)));
}

// void Aladdin_Ilator::BeginValidDecodeFuncDef(Aladdin_Ilator::CFunc* func,
// StrBuff& buff) const {
//   ILA_ASSERT(func->args.empty()); // no definition for uninterpreted funcs

//   StrBuff paramBuff;
//   AddParametersToBuffer(paramBuff, true);

//   auto type = (func->ret) ? GetCType(func->ret) : GetCType(func->ret_type);
//   auto args = (func->target)
//                   ? fmt::format("{} tmp_memory", GetCType(func->target))
//                   : "";

//   fmt::format_to(buff, "{return_type} {func_name}({argument}) {{\n",
//                  fmt::arg("return_type", type),
//                  fmt::arg("project", GetProjectName()),
//                  fmt::arg("func_name", func->name),
//                  fmt::arg("argument", fmt::to_string(paramBuff)));
// }

void Aladdin_Ilator::BeginMacroDef(Aladdin_Ilator::CFunc* func,
                                   StrBuff& buff) const {
  ILA_ASSERT(func->args.empty()); // no definition for uninterpreted funcs

  auto type = (func->ret) ? "output" : "";

  fmt::format_to(buff, "#define {func_name}({args}) {{\\\n",
                 fmt::arg("project", GetProjectName()),
                 fmt::arg("func_name", func->name), fmt::arg("args", type));
}

void Aladdin_Ilator::EndFuncDef(Aladdin_Ilator::CFunc* func,
                                StrBuff& buff) const {
  if (func->ret) {
    fmt::format_to(buff, "return {};\n\n", GetCName(func->ret));
  }
  fmt::format_to(buff, "}}\n\n");
}

void Aladdin_Ilator::EndMacroDef(Aladdin_Ilator::CFunc* func,
                                 StrBuff& buff) const {
  if (func->ret) {
    fmt::format_to(buff, "output = {};\\\n", GetCName(func->ret));
  }
  fmt::format_to(buff, "}}\n");
}

void Aladdin_Ilator::WriteFuncDecl(Aladdin_Ilator::CFunc* func,
                                   StrBuff& buff) const {
  auto type = (func->ret) ? GetCType(func->ret) : GetCType(func->ret_type);
  auto args = (func->target)
                  ? fmt::format("{}& tmp_memory", GetCType(func->target))
                  : "";

  if (!func->args.empty()) { // uninterpreted func only
    ILA_NOT_NULL(func->ret_type);
    std::vector<std::string> arg_list;
    for (const auto& a : func->args) {
      arg_list.push_back(GetCType(a));
    }
    args = fmt::format("{}", fmt::join(arg_list, ", "));
  } else {
    StrBuff paramBuff;
    AddParametersToBuffer(paramBuff, true, ",\n  ");
    args = fmt::to_string(paramBuff);
  }

  fmt::format_to(buff, "{return_type} {func_name}({argument});\n\n",
                 fmt::arg("return_type", type),
                 fmt::arg("func_name", func->name), fmt::arg("argument", args));
}

void Aladdin_Ilator::CommitSource(const std::string& file_name,
                                  const std::string& dir, const StrBuff& buff) {
  auto file_path = os_portable_append_dir(dir, file_name);

  if (file_path.find(".c") != std::string::npos) {
    sourcesList += " " + file_path.substr(file_path.rfind("/") + 1);
  }

  WriteFile(file_path, buff);
}

std::string Aladdin_Ilator::GetCxxType(const SortPtr& sort) {
  if (!sort) {
    return "void";
  } else if (sort->is_bool()) {
    return "bool";
  } else if (sort->is_bv()) {
    return fmt::format("sc_biguint<{}>", sort->bit_width());
  } else {
    ILA_ASSERT(sort->is_mem());
#ifdef ALADDIN_ILATOR_PRECISE_MEM
    return fmt::format(
        "std::map<sc_biguint<{addr_width}>, sc_biguint<{data_width}>>",
        fmt::arg("addr_width", sort->addr_width()),
        fmt::arg("data_width", sort->data_width()));
#else
    return "std::unordered_map<int, int>";
#endif
  }
}

std::string Aladdin_Ilator::GetCxxName(const ExprPtr& expr) {
  if (expr->is_var()) {
    return fmt::format("{}_{}", expr->host()->name().str(), expr->name().str());
  } else {
    return fmt::format("univ_var_{}", GetPivotalId(expr->name().id()));
  }
}

std::string Aladdin_Ilator::GetCType(const SortPtr& sort) {
  if (!sort) {
    return "void";
  } else if (sort->is_bool()) {
    return "bool";
  } else if (sort->is_bv()) {
    if (sort->bit_width() > 1) {
      return fmt::format("_BitInt({})", sort->bit_width());
    } else {
      return "bool";
    }
  } else {
    ILA_ASSERT(sort->is_mem());
    uint64_t num_elements = 1 << sort->addr_width();

    if (sort->data_width() > 1) {
      return fmt::format("_BitInt({data_width})[{num_elements}]",
                         fmt::arg("num_elements", num_elements),
                         fmt::arg("data_width", sort->data_width()));
    } else {
      return fmt::format("bool[{num_elements}]",
                         fmt::arg("num_elements", num_elements));
    }
  }
}

std::string Aladdin_Ilator::GetCName(const ExprPtr& expr) {
  if (expr->is_var()) {
    return fmt::format("{}_{}", expr->host()->name().str(), expr->name().str());
  } else {
    return fmt::format("univ_var_{}", GetPivotalId(expr->name().id()));
  }
}

std::string Aladdin_Ilator::GetCNameWithType(const ExprPtr& expr) {

  std::string name = "";
  if (expr->is_var()) {
    name = fmt::format("{}_{}", expr->host()->name().str(), expr->name().str());
  } else {
    name = fmt::format("univ_var_{}", GetPivotalId(expr->name().id()));
  }

  const SortPtr sort = expr->sort();

  if (!sort) {
    return "void";
  } else if (sort->is_bool()) {
    return "bool " + name;
  } else if (sort->is_bv()) {
    if (sort->bit_width() > 1) {
      return fmt::format("_BitInt({}) {}", sort->bit_width(), name);
    } else {
      return fmt::format("bool {}", name);
    }
  } else {
    ILA_ASSERT(sort->is_mem());
    uint64_t num_elements = (uint64_t)1 << sort->addr_width();

    if (sort->data_width() > 1) {
      return fmt::format("_BitInt({data_width}) {name}[{num_elements}]",
                         fmt::arg("num_elements", num_elements),
                         fmt::arg("name", name),
                         fmt::arg("data_width", sort->data_width()));
    } else {
      return fmt::format("bool {name}[{num_elements}]",
                         fmt::arg("num_elements", num_elements),
                         fmt::arg("name", name));
    }
  }
}

uint64_t Aladdin_Ilator::GetWordSize(const ExprPtr& expr) {

  const SortPtr sort = expr->sort();

  if (sort->is_bool()) {
    return 1;
  }
  if (sort->is_bv()) {
    return ((sort->bit_width() - 1) / 8) + 1;
  }
  return ((sort->data_width() - 1) / 8) + 1;
}

uint64_t Aladdin_Ilator::GetNumBytes(const ExprPtr& expr) {

  const SortPtr sort = expr->sort();

  if (sort->is_bool() || sort->is_bv()) {
    return GetWordSize(expr);
  }

  std::cout << ((uint64_t)1 << ((uint64_t)sort->addr_width())) << '\n';

  // Frankly if it is expecting to use more than 2^64 bytes it will have other
  // problems
  return (((uint64_t)1 << ((uint64_t)sort->addr_width())) * GetWordSize(expr));
}

std::string Aladdin_Ilator::GetValidFuncName(const InstrLvlAbsCnstPtr& m) {
  return fmt::format("valid_{host}", fmt::arg("host", m->name().str()));
}

std::string Aladdin_Ilator::GetDecodeFuncName(const InstrPtr& instr) {
  return fmt::format("decode_{host}_{instr}",
                     fmt::arg("host", instr->host()->name().str()),
                     fmt::arg("instr", instr->name().str()));
}

std::string Aladdin_Ilator::GetUpdateFuncName(const InstrPtr& instr) {
  return fmt::format("update_{host}_{instr}",
                     fmt::arg("host", instr->host()->name().str()),
                     fmt::arg("instr", instr->name().str()));
}

std::string Aladdin_Ilator::GetMemoryFuncName(const ExprPtr& expr) {
  ILA_ASSERT(expr->is_mem());
  if (asthub::GetUidExprOp(expr) == AstUidExprOp::kIfThenElse) {
    return fmt::format("ite_{}", GetPivotalId(expr->name().id()));
  } else {
    return fmt::format("store_{}", GetPivotalId(expr->name().id()));
  }
}

} // namespace ilang
