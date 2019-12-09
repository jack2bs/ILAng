/// \file CounterExample Extractor
// ---Hongce Zhang

#include <ilang/config.h>

#include <ilang/util/container_shortcut.h>
#include <ilang/util/log.h>
#include <ilang/util/str_util.h>
#include <ilang/vtarget-out/inv-syn/cex_extract.h>

#include <vcdparser/VCDFileParser.hpp>

#include <fstream>
#include <sstream>

namespace ilang {

static std::string val2str(const VCDValue& v) {
  std::stringstream ret;

  switch (v.get_type()) {
  case (VCD_SCALAR):
    ret << "1'b" << VCDValue::VCDBit2Char(v.get_value_bit());
    break;
  case (VCD_VECTOR): {
    const VCDBitVector* vecval = v.get_value_vector();
    ret << std::to_string(vecval->size()) << "'b";
    for (auto it = vecval->begin(); it != vecval->end(); ++it)
      ret << VCDValue::VCDBit2Char(*it);
  } break;
  case (VCD_REAL):
    ret << v.get_value_real();
    break;
  default:
    ILA_ERROR << "Unknown value type!";
  }
  return ret.str();
}

std::string prepend(const std::string& prefix, const std::string& sig_name) {
  if (prefix.empty())
    return sig_name;
  return prefix + "." + sig_name;
}

/*
bool in_scope(VCDScope * sc, const std::string & sc_name) {
  VCDScope * sc_traverse = sc;
  while(sc_traverse) {
    if (sc_traverse->name == sc_name)
      return true;
    sc_traverse = sc_traverse->parent;
  }
  return false;
}
*/

std::string collect_scope(VCDScope* sc) {
  std::string ret;
  while (sc) {
    ret = sc->name + "." + ret;
    sc = sc->parent;
  }
  return ret;
}

void CexExtractor::parse_from(const std::string& vcd_file_name,
                              const std::string& scope, is_reg_t is_reg,
                              bool reg_only) {

  cex.clear();

  VCDFileParser parser;
  VCDFile* trace = parser.parse_file(vcd_file_name);

  if (!trace) {
    ILA_ERROR << "Error while reading waveform from: " << vcd_file_name;
    return;
  }

  ILA_NOT_NULL(trace->get_scope("$root"));
  VCDScope* top = NULL;
  for (VCDScope* c : trace->get_scope("$root")->children)
    if (c->name == "top")
      top = c;
  ILA_NOT_NULL(top);

  // determine the start signal time
  std::string start_sig_hash;
  for (VCDSignal* root_sig : top->signals) {
    // find the hash of it
    if (root_sig->reference == "__START__" ||
        root_sig->reference == "__START__[0:0]")
      start_sig_hash = root_sig->hash;
  }
  if (start_sig_hash.empty()) {
    ILA_ERROR << "Error analyzing waveform. "
              << "It is not a trace generated by Verilog Verification Target "
                 "Generation. "
              << "__START__ signal not found in $root scope";
    return;
  }

  VCDSignalValues* start_sig_vals = trace->get_signal_value(start_sig_hash);
  VCDTime start_time = -1;
  for (VCDTimedValue* tv : *start_sig_vals) {
    if (val2str(*(tv->value)) == "1'b1") {
      start_time = tv->time;
      break;
    }
  }

  if (start_time == -1) {
    ILA_ERROR << "Start time not found from waveform!";
    return;
  }

  std::vector<VCDSignal*>* sigs = trace->get_signals();

  for (VCDSignal* sig : *sigs) {

    // ensure it is only register
    if (sig->type != VCDVarType::VCD_VAR_REG)
      continue;

    // check scope
    // if (! in_scope(sig->scope, scope))
    //  continue;

    auto scopes = collect_scope(sig->scope);

    // check scope -- only the top level
    if (!(StrStartsWith(scopes, "$root.top." + scope + ".") ||
          StrStartsWith(scopes, scope + ".")))
      continue;

    auto vlg_name = ReplaceAll(scopes + sig->reference, "$root.top.", "");

    std::string check_name = vlg_name;
    {
      auto pos = check_name.find('[');
      if (pos != std::string::npos)
        check_name = check_name.substr(0, pos);
    }

    bool is_this_var_reg = is_reg(check_name);

    if (reg_only && !is_this_var_reg)
      continue;

    auto vlg_val_ptr = trace->get_signal_value_at(sig->hash, start_time);

    if (vlg_val_ptr == nullptr) {
      ILA_WARN << "Parsing VCD: " << vlg_name << " gets Xs. Ignored.";
      continue;
    }

    std::string val = val2str(*vlg_val_ptr);

    cex.insert(std::make_pair(vlg_name, val));
    cex_is_reg.insert(std::make_pair(vlg_name, is_this_var_reg));

  } // for sig

  ILA_ASSERT(!cex.empty()) << "No counterexample is extracted!";

} // parse_from

CexExtractor::CexExtractor(const std::string& vcd_file_name,
                           const std::string& scope, is_reg_t is_reg,
                           bool reg_only) {
  parse_from(vcd_file_name, scope, is_reg, reg_only);
}

/// create from a existing file
CexExtractor::CexExtractor(const std::string& fn) {
  std::ifstream fin(fn);
  if (!fin.is_open()) {
    ILA_ERROR << "Failed to read : " << fn;
    return;
  }
  unsigned pairs;
  fin >> pairs;

  std::string name, val;
  for (unsigned idx = 0; idx < pairs; ++idx) {
    fin >> name;
    fin >> val;
    if (name.empty() || val.empty()) {
      ILA_ERROR << "Failed to read cex from file!";
      break;
    }
    cex.insert(std::make_pair(name, val));
    cex_is_reg.insert(std::make_pair(name, true));
  }
}

// -------------------- MEMBERS ------------------ //
/// return a string to be added to the design
std::string
CexExtractor::GenInvAssert(const std::string& prefix,
                           const std::set<std::string>& focus_name) const {

  std::string ret = "(1'b1 == 1'b1)"; // true
  for (auto&& nv : cex) {
    if (!cex_is_reg.at(nv.first))
      continue;
    auto fullname = prepend(prefix, ReplaceAll(nv.first, "[0:0]", ""));
    auto check_name = fullname; // remove further of [][:]
    auto pos = check_name.rfind('[');
    if (pos != std::string::npos)
      check_name = check_name.substr(0, pos);
    if (!focus_name.empty() && !IN(check_name, focus_name))
      continue;
    ret += "\n&& (" + fullname + " == " + nv.second + ")";
  }
  return ret;
}

const CexExtractor::cex_t& CexExtractor::GetCex() const { return cex; }

// save to file
void CexExtractor::StoreCexToFile(const std::string& fn, const cex_t& c) {
  std::ofstream fout(fn);

  if (!fout.is_open()) {
    ILA_ERROR << "Failed to read : " << fn;
    return;
  }

  fout << c.size() << "\n";
  for (auto&& nv : c) {
    ILA_ERROR_IF(S_IN(' ', nv.first) || S_IN('\r', nv.first) ||
                 S_IN('\t', nv.first) || S_IN('\n', nv.first))
        << nv.first << " contains space!";
    fout << nv.first << " " << nv.second << std::endl;
  }
}
// save to file (invoke within)
void CexExtractor::StoreCexToFile(const std::string& fn) const {
  StoreCexToFile(fn, cex);
}

// generalize cex
void CexExtractor::DropStates(const std::vector<std::string>& vnames) {
  for (auto&& n : vnames) {
    if (IN(n, cex)) {
      cex.erase(n);
      cex_is_reg.erase(n);
    }
  }
}

}; // namespace ilang
