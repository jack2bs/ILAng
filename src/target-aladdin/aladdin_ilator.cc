/// \file
/// Implementation of the class Aladdin_Ilator.

#include <ilang/target-aladdin/aladdin_ilator.h>

#include <fstream>
#include <math.h>

#include <fmt/format.h>
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

static std::string sourcesList = "";
static void WriteFile(const std::string& file_path, const fmt::memory_buffer& buff) {
  if (file_path.find(".c") != std::string::npos)
  {
    sourcesList += " " + file_path.substr(file_path.rfind("/") + 1);
  }
  
  
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

  // GET ALL OF THE MEM VARIABLES TYPES (CACHE/DMA/SCRATCHPAD)
  size_t numStates = m_->state_num();
  for (size_t i = 0; i < numStates; i++)
  {
    auto s = m_->state(i);
    if (s->sort()->is_mem())
    {
      memory_types.insert({s->name().str(), spad});
    }
  }
  
  // size_t numStates = m_->state_num();

  // for (size_t i = 0; i < numStates; i++)
  // {
  //   auto s = m_->state(i);
  //   if (s->sort()->is_mem())
  //   {
  //     std::cout << s->name();
  //   }
  // }
  
}

Aladdin_Ilator::~Aladdin_Ilator() { Reset(); }

void Aladdin_Ilator::Generate(const std::string& dst, bool opt) {
  // sanity checks and initialize
  if (!SanityCheck() || !Bootstrap(dst, opt)) {
    return;
  }

  auto status = true;
  ILA_INFO << "Start generating SystemC simulator of " << m_;

  // non-instruction basics
  status &= GenerateIlaBasics(os_portable_append_dir(dst, kDirSrc));

  // instruction semantics (decode and updates)
  for (auto& instr : absknob::GetInstrTree(m_)) {
    status &= GenerateInstrContent(instr, os_portable_append_dir(dst, kDirSrc));
  }

  // memory updates
  status &= GenerateMemoryUpdate(os_portable_append_dir(dst, kDirSrc));

  // constant memory
  status &= GenerateConstantMemory(os_portable_append_dir(dst, kDirSrc));

  // initial condition setup
  status &= GenerateInitialSetup(os_portable_append_dir(dst, kDirSrc));

  // execution kernel
  status &= GenerateExecuteKernel(os_portable_append_dir(dst, kDirSrc));

  // shared header (input, state, func., etc.)
  status &= GenerateGlobalHeader(os_portable_append_dir(dst, kDirSrc));

  // file to define the inputs
  status &= GenerateInputTemplate(os_portable_append_dir(dst, kDirSrc));

  // file which configures the variables
  status &= GenerateConfigTemplate(os_portable_append_dir(dst, kDirSrc));

  // cmake support, e.g., recipe and templates
  status &= GenerateBuildSupport(os_portable_append_dir(dst, kDirSrc));




  // clean up if something went wrong
  if (status) {
    ILA_INFO << "Sucessfully generate SystemC simulator at " << dst;
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

  // memory updates
  for (auto f : memory_updates_) {
    delete f.second;
  }
  memory_updates_.clear();

  // others
  source_files_.clear();
  const_mems_.clear();
  global_vars_.clear();
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
    auto valid_func = RegisterFunction(GetValidFuncName(m), valid_expr);
    BeginFuncDef(valid_func, buff);
    ExprVarMap lut;
    ILA_CHECK(RenderExpr(valid_expr, buff, lut));
    fmt::format_to(buff, "{var} = {local_name};\n",
                   fmt::arg("var", GetCNameWithType(valid_expr)),
                   fmt::arg("local_name", LookUp(valid_expr, lut)));
    EndFuncDef(valid_func, buff);
  };

  // traverse the hierarchy
  m_->DepthFirstVisit(PerIla);

  // record and write to file
  CommitSource("all_valid_funcs_in_hier.c", dir, buff);
  return true;
}

bool Aladdin_Ilator::GenerateInstrContent(const InstrPtr& instr,
                                  const std::string& dir) {
  StrBuff buff;
  ExprVarMap lut;

  // include headers
  fmt::format_to(buff, "#include \"{}.h\"\n", GetProjectName());

  // decode function
  auto decode_expr = instr->decode();
  auto decode_func = RegisterFunction(GetDecodeFuncName(instr), decode_expr);
  BeginFuncDef(decode_func, buff);
  lut.clear();
  if (!RenderExpr(decode_expr, buff, lut)) {
    return false;
  }
  fmt::format_to(buff, "{var} = {local_name};\n",
                 fmt::arg("var", GetCNameWithType(decode_expr)),
                 fmt::arg("local_name", LookUp(decode_expr, lut)));
  EndFuncDef(decode_func, buff);

  // next state
  auto update_func = RegisterFunction(GetUpdateFuncName(instr));
  BeginFuncDef(update_func, buff);
  lut.clear();
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
      if (!RenderExpr(update_expr, buff, lut)) {
        return false;
      }
      fmt::format_to(buff, "{type_name} {local_var}_nxt_holder = {local_var};\n",
                     fmt::arg("type_name", GetCType(update_expr)),
                     fmt::arg("local_var", LookUp(update_expr, lut)));
    } else { // memory (one copy for performance, require special handling)
      if (HasLoadFromStore(update_expr)) {
        // Load from store is not allowed
        return false;
      }
      if (!RenderExpr(update_expr, buff, lut)) {
        return false;
      }
      fmt::format_to(buff,
                     "_BitInt({mem_addr_width}) {local_var}_addr_nxt_holder = {local_var}_addr;\n"
                     "_BitInt({mem_data_width}) {local_var}_data_nxt_holder = {local_var}_data;\n",
                     fmt::arg("mem_addr_width", update_expr->sort()->addr_width()),
                     fmt::arg("mem_data_width", update_expr->sort()->data_width()),
                     fmt::arg("local_var",  LookUp(update_expr, lut)));
    }
  }
  // update state
  for (auto& s : updated_states) {
    auto curr = instr->host()->state(s);
    auto next = instr->update(s);
    if (!curr->is_mem()) {
      fmt::format_to(buff, "{current} = {next_value}_nxt_holder;\n",
                     fmt::arg("current", GetCName(curr)),
                     fmt::arg("next_value", LookUp(next, lut)));
    } else {
      // UPDATE ANY CHANGED VALUES
      fmt::format_to(buff,
                     "{current}[{next_value}_addr_nxt_holder] = {next_value}_data_nxt_holder;\n",
                     fmt::arg("current", GetCName(curr)),
                     fmt::arg("next_value", LookUp(next, lut)));
    }
  }

  // // add update states logging
  // fmt::format_to(buff, "#ifdef Aladdin_Ilator_VERBOSE\n");
  // fmt::format_to(buff, "instr_update_log << \"No.\" << std::dec << GetInstrCntr() << '\\t' << "
  //                 "\"{instr_name} state updates:\" << std::endl;\n",
  //                 fmt::arg("instr_name", instr->name().str()));
  // for (auto& s : updated_states) {
  //   auto curr = instr->host()->state(s);
  //   if (!curr->is_mem()) {
  //     fmt::format_to(buff,
  //                    "instr_update_log << \"    {state_name} => \" << "
  //                    "std::hex << \"0x\" << {state_name} << std::endl; \n",
  //                    fmt::arg("state_name", GetCName(curr)));
  //   } else {
  //     fmt::format_to(buff,
  //                    "instr_update_log << \"    {state_name} get updated\" "
  //                    "<< std::endl;\n",
  //                    fmt::arg("state_name", GetCName(curr)));
  //   }
  // }
  // fmt::format_to(buff, "instr_update_log << std::endl;\n");
  // fmt::format_to(buff, "#endif\n");

  EndFuncDef(update_func, buff);

  // record and write to file
  CommitSource(fmt::format("idu_{}.c", instr->name().str()), dir, buff);
  return true;
}

bool Aladdin_Ilator::GenerateMemoryUpdate(const std::string& dir) {

  // helper for traversing memory updates
  class MemUpdateVisiter {
  public:
    MemUpdateVisiter(Aladdin_Ilator* h, StrBuff& b, ExprVarMap& l)
        : host(h), buff_ref(b), lut_ref(l) {}

    bool pre(const ExprPtr& expr) {
      // stop traversing when reaching memory ITE (stand-alone func)
      if (expr->is_mem() && expr->is_op() &&
          asthub::GetUidExprOp(expr) == AstUidExprOp::kIfThenElse) {
        host->DfsExpr(expr, buff_ref, lut_ref);
        return true;
      } else {
        return false;
      }
    }

    void post(const ExprPtr& expr) { host->DfsExpr(expr, buff_ref, lut_ref); }

    Aladdin_Ilator* host;
    StrBuff& buff_ref;
    ExprVarMap lut_ref;
  };

  auto RenderMemUpdate = [this](const ExprPtr& e, StrBuff& b, ExprVarMap& l) {
    auto mem_visiter = MemUpdateVisiter(this, b, l);
    e->DepthFirstVisitPrePost(mem_visiter);
  };

  // helpers for managing files
  int file_cnt = 0;
  auto GetMemUpdateFile = [&file_cnt]() {
    return fmt::format("memory_update_functions_{}.c", file_cnt++);
  };

  StrBuff buff;
  auto StartNewFile = [this, &buff]() {
    buff.clear();
    fmt::format_to(buff, "#include \"{}.h\"\n", GetProjectName());
  };

  // start generating
  StartNewFile();
  ExprVarMap lut;

  for (auto mem_update_func_it : memory_updates_) {
    auto& mem_update_func = mem_update_func_it.second;
    ILA_NOT_NULL(mem_update_func);
    auto& mem = mem_update_func->target;

    lut.clear();

    BeginFuncDef(mem_update_func, buff);

    if (asthub::GetUidExprOp(mem) == AstUidExprOp::kStore) {
      RenderMemUpdate(mem, buff, lut);
    } else { // ite
      RenderExpr(mem->arg(0), buff, lut);
      auto lut_local_true = lut;
      auto& lut_local_false = lut; // reuse

      fmt::format_to(buff, "if ({}) {{\n", LookUp(mem->arg(0), lut));
      RenderMemUpdate(mem->arg(1), buff, lut_local_true);
      fmt::format_to(buff, "}} else {{\n");
      RenderMemUpdate(mem->arg(2), buff, lut_local_false);
      fmt::format_to(buff, "}}\n");
    }

    EndFuncDef(mem_update_func, buff);

    if (buff.size() > 50000) {
      CommitSource(GetMemUpdateFile(), dir, buff);
      StartNewFile();
    }
  }

  CommitSource(GetMemUpdateFile(), dir, buff);
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
      addr_data_pairs.push_back(fmt::format("{mem_name}[{addr}] = {data};",
                                            fmt::arg("mem_name", GetCName(mem)),
                                            fmt::arg("addr", it.first),
                                            fmt::arg("data", it.second)));
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
  auto init_func = RegisterFunction("setup_initial_condition");
  StrBuff buff;
  fmt::format_to(buff, "#include \"{}.h\"\n", GetProjectName());
  BeginFuncDef(init_func, buff);
  for (auto pair : init_values) {
    fmt::format_to(buff, "{var_name} = {var_value};\n",
                   fmt::arg("var_name", GetCName(pair.first)),
                   fmt::arg("var_value", pair.second));
  }
  EndFuncDef(init_func, buff);
  CommitSource("setup_initial_condition.c", dir, buff);
  return true;
}

bool Aladdin_Ilator::GenerateExecuteKernel(const std::string& dir) {
  StrBuff buff;

  fmt::format_to(buff, "#include \"{}.h\"\n", GetProjectName());
  fmt::format_to(buff, "#include \"gem5/dma_interface.h\"\n");
  



  auto kernel_func = RegisterFunction("compute");
  BeginFuncDef(kernel_func, buff);

  // setup initial condition
  fmt::format_to(buff, 
                       "  int instr_cntr = 0;\n"
                       "  setup_initial_condition();\n\n"
                       "  main_while: while (instr_cntr < cycles_to_simulate) {{\n\n");

  // read in input value
  for (size_t i = 0; i < m_->input_num(); i++) {
    fmt::format_to(buff, "    dmaLoad(&{input_name}, {input_name}_inps + instr_cntr, sizeof({input_name}));\n",
                   fmt::arg("input_name", GetCName(m_->input(i))));
  }

  // instruction execution
  auto ExecInstr = [this, &buff](const InstrPtr& instr, bool child) {
    fmt::format_to(
        buff,
      "\n{ind}    if ({valid_func_name}() && {decode_func_name}()) {{\n"
        "{ind}      {update_func_name}();\n"
        "{ind}{child_counter}"
        "{ind}    }}\n\n",
        fmt::arg("ind", (child ? "  " : "")),
        fmt::arg("valid_func_name", GetValidFuncName(instr->host())),
        fmt::arg("decode_func_name", GetDecodeFuncName(instr)),
        fmt::arg("update_func_name", GetUpdateFuncName(instr)),
        fmt::arg("child_counter", (child ? "      schedule_counter++;\n" : "")),
        fmt::arg("instr_name", instr->name().str()));
  };

  auto top_instrs = absknob::GetInstr(m_);
  auto all_instrs = absknob::GetInstrTree(m_);

  // top-level instr
  for (auto& instr : top_instrs) {
    ExecInstr(instr, false);
  }

  // child instr
  fmt::format_to(buff, "    while (1) {{\n"
                       "      int schedule_counter = 0;\n");
  std::set<InstrPtr> tops(top_instrs.begin(), top_instrs.end());
  for (auto& instr : all_instrs) {
    if (tops.find(instr) == tops.end()) {
      ExecInstr(instr, true);
    }
  }
  fmt::format_to(buff, "      if (schedule_counter == 0) {{\n"
                       "        break;\n"
                       "      }}\n"
                       "    }}\n"
                       "    instr_cntr++;\n"
                       "  }}\n");

  // done
  EndFuncDef(kernel_func, buff);

  CommitSource("compute.c", dir, buff);
  return true;
}

bool Aladdin_Ilator::GenerateGlobalHeader(const std::string& dir) {
  StrBuff buff;

  fmt::format_to(buff,
                 "#ifndef _{project}_H_\n"
                 "#define _{project}_H_\n",
                 fmt::arg("project", GetProjectName()));

  fmt::format_to(buff, "\n#include <stdbool.h>\n");

  // function declarations
  fmt::format_to(buff, "\n// function declarations\n");
  for (auto& func : functions_) {
    WriteFuncDecl(func.second, buff);
  }
  for (auto& func : externs_) {
    WriteFuncDecl(func.second, buff);
  }
  for (auto& func : memory_updates_) {
    WriteFuncDecl(func.second, buff);
  }

  // inputs
  fmt::format_to(buff, "\n// inputs\n");
  for (auto& var : absknob::GetInp(m_)) {
    fmt::format_to(buff,
                   "extern {type} {name};\n"
                   "extern {type} * {name}_inps;\n",
                   fmt::arg("type", GetCType(var)),
                   fmt::arg("name", GetCName(var)));
  }

  // state and global vars
  fmt::format_to(buff, "\n// states\n");
  for (auto& var : absknob::GetSttTree(m_)) {
    fmt::format_to(buff,
                   "extern {var};\n",
                   fmt::arg("var", GetCNameWithType(var)));
  }
  fmt::format_to(buff, "\n// global vars\n");
  for (auto& var : global_vars_) {
    fmt::format_to(buff,
                   "extern {var};\n",
                   fmt::arg("var", GetCNameWithType(var)));
  }
  fmt::format_to(buff,
                 "extern const int cycles_to_simulate;\n");

  // memory constant
  fmt::format_to(buff, "\n// memory constants\n");
  for (auto& mem : const_mems_) {
    fmt::format_to(buff, 
                  "extern static {var};\n",
                   fmt::arg("var", GetCNameWithType(mem)));
  }

  fmt::format_to(buff,
                 "\n#endif//_{project}_H_\n",
                 fmt::arg("project", GetProjectName()));

  // write to file
  auto file_path = os_portable_append_dir(dir, GetProjectName() + ".h");
  WriteFile(file_path, buff);
  return true;
}

bool Aladdin_Ilator::GenerateBuildSupport(const std::string& dir) {
  
  StrBuff buff;
  buff.clear();

  fmt::format_to(buff, "#include \"{}.h\"\n", GetProjectName());

  // input
  if (absknob::GetInp(m_).size())
  {
    fmt::format_to(buff, "\n// inputs\n");
  }
  for (auto& var : absknob::GetInp(m_)) {
    fmt::format_to(buff,
                   "{var};\n",
                   fmt::arg("var", GetCNameWithType(var)));
  }

  // state and global vars
  if (absknob::GetSttTree(m_).size())
  {
    fmt::format_to(buff, "\n// states\n");
  }
  for (auto& var : absknob::GetSttTree(m_)) {
    fmt::format_to(buff, "{var};\n",
                   fmt::arg("var", GetCNameWithType(var)));
  }

  if (global_vars_.size())
  {  
    fmt::format_to(buff, "\n// global vars\n"); 
  }
  for (auto& var : global_vars_) {
    fmt::format_to(buff, "{var};\n",
                   fmt::arg("var", GetCNameWithType(var)));
  }

  // memory constant
  if (const_mems_.size())
  {
    fmt::format_to(buff, "\n// memory constants\n");
  }
  for (auto& mem : const_mems_) {
    fmt::format_to(buff, "static {var};\n",
                   fmt::arg("var", GetCNameWithType(mem)));
  }

  // main function
  static const char* kSimEntryTemplate =
      "\nint main(int argc, char* argv[]) {{\n"
      "  compute();\n"
      "  return 0; \n"
      "}}\n";

  fmt::format_to(buff, kSimEntryTemplate);
      
  WriteFile(os_portable_append_dir(dir, "main.c"), buff);


  buff.clear();

  fmt::format_to(buff, 
                "SRCS={src_text}\n\n"
                "ACCEL_NAME = ila\n"
                "TEST_BIN = $(ACCEL_NAME)\n"
                "export TRACE_OUTPUT_DIR=$(ACCEL_NAME)\n"
                "ifndef WORKLOAD\n"
                "  export WORKLOAD=compute\n"
                "endif\n"
                "include ../common/Makefile.common\n"
                "include ../common/Makefile.tracer\n",
                fmt::arg("src_text", sourcesList));

  WriteFile(os_portable_append_dir(dir, "Makefile"), buff);

  return true;
}

bool Aladdin_Ilator::GenerateInputTemplate(const std::string& dir) {

  StrBuff buff;

  fmt::format_to(buff, "#include \"{}.h\"\n\n", GetProjectName());


  fmt::format_to(buff, "// How many cycles would you like to simulate?\n"
                       "#define CYCLES_TO_SIMULATE 5\n\n"
                       "const int cycles_to_simulate = CYCLES_TO_SIMULATE;\n\n");

  for (auto& var : absknob::GetInp(m_)) {
    fmt::format_to(buff,
                   "// Please enter values for the input \"{var}\"\n"
                   "{c_name}_inps_vals[CYCLES_TO_SIMULATE] = {{\n"
                   "  {def1},\n"
                   "  {def2},\n"
                   "  {def1},\n"
                   "  {def2},\n"
                   "  {def1},\n"
                   "}};\n\n"
                   "{type} * {name}_inps = &{name}_inps_vals[0];\n\n",
                   fmt::arg("var", var->name().str()),
                   fmt::arg("c_name", GetCNameWithType(var)),
                   fmt::arg("def1", var->is_bool() ? "false" : "0"),
                   fmt::arg("def2", var->is_bool() ? "true" : "1"),
                   fmt::arg("type", GetCType(var)),
                   fmt::arg("name", GetCName(var)));
  }

  auto entry_path =
      os_portable_append_dir(dir, "input_vals.c");
  
  WriteFile(entry_path, buff);
  return true;

}

void Aladdin_Ilator::AddConfigLineToBuff(const ilang::ExprPtr & var, StrBuff & buff) {
  if (!var->is_mem())
  {    
    fmt::format_to(buff,
                 "partition,complete,{var},{num_bytes},{word_size}\n",
                 fmt::arg("var", GetCName(var)),
                 fmt::arg("num_bytes", GetNumBytes(var)),
                 fmt::arg("word_size", GetWordSize(var)));
  }
  else
  {
    switch (memory_types.at(var->name().str())) {
    case spad:
      fmt::format_to(buff,
                 "partition,cyclic,{var},{num_bytes},{word_size},1\n",
                 fmt::arg("var", GetCName(var)),
                 fmt::arg("num_bytes", GetNumBytes(var)),
                 fmt::arg("word_size", GetWordSize(var)));
      break;
    case dma:
      break;
    case cache:
      fmt::format_to(buff,
                 "cache,{var},{word_size}\n",
                 fmt::arg("var", GetCName(var)),
                 fmt::arg("word_size", GetWordSize(var)));
      break;
    }
  }
} 

bool Aladdin_Ilator::GenerateConfigTemplate(const std::string& dir) {
  
  StrBuff buff;

  // inputs
  // fmt::format_to(buff, "\n// inputs\n");
  for (auto& var : absknob::GetInp(m_)) {
    AddConfigLineToBuff(var, buff);
    fmt::format_to(buff,
                 "partition,complete,{var}_inps,8,8\n",
                 fmt::arg("var", GetCName(var)),
                 fmt::arg("num_bytes", GetNumBytes(var)),
                 fmt::arg("word_size", GetWordSize(var)));
  }

  // state and global vars
  // fmt::format_to(buff, "\n// states\n");
  for (auto& var : absknob::GetSttTree(m_)) {
    AddConfigLineToBuff(var, buff);
  }
  // fmt::format_to(buff, "\n// global vars\n");
  for (auto& var : global_vars_) {
    AddConfigLineToBuff(var, buff);
  }
  fmt::format_to(buff,
                 "partition,complete,cycles_to_simulate,4,4\n");

  // memory constant
  // fmt::format_to(buff, "\n// memory constants\n");
  for (auto& mem : const_mems_) {
    AddConfigLineToBuff(mem, buff);
  }


  fmt::format_to(buff,
                  "cycle_time,6\n"
                  "unrolling,compute,main_while,1\n");


  auto entry_path =
      os_portable_append_dir(dir, "config");
  
  WriteFile(entry_path, buff);
  return true;
}

bool Aladdin_Ilator::RenderExpr(const ExprPtr& expr, StrBuff& buff, ExprVarMap& lut) {

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

Aladdin_Ilator::CFunc* Aladdin_Ilator::RegisterFunction(const std::string& func_name,
                                          ExprPtr return_expr) {
  auto func = new CFunc(func_name, return_expr);
  auto [it, status] = functions_.insert({func->name, func});
  ILA_ASSERT(status);
  return func;
}

Aladdin_Ilator::CFunc* Aladdin_Ilator::RegisterExternalFunc(const FuncPtr& func) {
  auto func_cxx = new CFunc(func->name().str(), func->out());
  auto [it, status] = externs_.insert({func_cxx->name, func_cxx});
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

Aladdin_Ilator::CFunc* Aladdin_Ilator::RegisterMemoryUpdate(const ExprPtr& mem) {
  auto func_cxx = new CFunc(GetMemoryFuncName(mem), NULL, mem);
  auto [it, status] = memory_updates_.insert({func_cxx->name, func_cxx});
  // memory updates can have multiple occurrence
  if (!status) {
    delete func_cxx;
  }
  return it->second;
}

void Aladdin_Ilator::BeginFuncDef(Aladdin_Ilator::CFunc* func, StrBuff& buff) const {
  ILA_ASSERT(func->args.empty()); // no definition for uninterpreted funcs

  auto type = (func->ret) ? GetCType(func->ret) : GetCType(func->ret_type);
  auto args = (func->target)
                  ? fmt::format("{} tmp_memory", GetCType(func->target))
                  : "";

  fmt::format_to(buff, "{return_type} {func_name}({argument}) {{\n",
                 fmt::arg("return_type", type),
                 fmt::arg("project", GetProjectName()),
                 fmt::arg("func_name", func->name),
                 fmt::arg("argument", args));
}

void Aladdin_Ilator::EndFuncDef(Aladdin_Ilator::CFunc* func, StrBuff& buff) const {
  if (func->ret) {
    fmt::format_to(buff, "return {};\n", GetCName(func->ret));
  }
  fmt::format_to(buff, "}}\n");
}

void Aladdin_Ilator::WriteFuncDecl(Aladdin_Ilator::CFunc* func, StrBuff& buff) const {
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
  }

  fmt::format_to(buff, "{return_type} {func_name}({argument});\n",
                 fmt::arg("return_type", type),
                 fmt::arg("func_name", func->name), fmt::arg("argument", args));
}

void Aladdin_Ilator::CommitSource(const std::string& file_name, const std::string& dir,
                          const StrBuff& buff) {
  auto file_path = os_portable_append_dir(dir, file_name);
  auto [it, ret] = source_files_.insert(file_name);
  ILA_ASSERT(ret) << "Duplicated source file name " << file_name;

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
    return fmt::format("_BitInt({})", sort->bit_width());
  } else {
    ILA_ASSERT(sort->is_mem());
    uint64_t num_elements = 1 << sort->addr_width();
    return fmt::format(
        "_BitInt({data_width})[{num_elements}]",
        fmt::arg("num_elements", num_elements),
        fmt::arg("data_width", sort->data_width()));
  }
}

std::string Aladdin_Ilator::GetCName(const ExprPtr& expr) {
  if (expr->is_var()) {
    return fmt::format("{}_{}", expr->host()->name().str(), expr->name().str());
  } else {
    return fmt::format("univ_var_{}", GetPivotalId(expr->name().id()));
  }
}

std::string Aladdin_Ilator::GetCNameWithType(const ExprPtr & expr) {
  
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
    return fmt::format("_BitInt({}) {}", sort->bit_width(), name);
  } else {
    ILA_ASSERT(sort->is_mem());
    uint64_t num_elements = 1 << sort->addr_width();
    return fmt::format(
        "_BitInt({data_width}) {name}[{num_elements}]",
        fmt::arg("num_elements", num_elements),
        fmt::arg("name", name),
        fmt::arg("data_width", sort->data_width()));
  }
}

uint64_t Aladdin_Ilator::GetWordSize(const ExprPtr & expr) {
  
  const SortPtr sort = expr->sort();

  if (sort->is_bool())
  {
    return 1;
  }
  if (sort->is_bv()) 
  {
    return ((sort->bit_width() - 1) / 8 ) + 1;
  }
  return ((sort->data_width() - 1) / 8 ) + 1;
}

uint64_t Aladdin_Ilator::GetNumBytes(const ExprPtr & expr) {
  
  const SortPtr sort = expr->sort();

  if (sort->is_bool() || sort->is_bv()) {
    return GetWordSize(expr);
  }

  // Frankly if it is expecting to use more than 2^64 bytes it will have other problems
  return ((1 << sort->addr_width()) * GetWordSize(expr));
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
