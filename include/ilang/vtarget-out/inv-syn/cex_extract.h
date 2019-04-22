/// \file The module to wrap vcd-parser (CounterExample Extractor)
// ---Hongce Zhang

#ifndef CEX_EXTRACT_H__
#define CEX_EXTRACT_H__

#include <string>
#include <map>
#include <functional>

namespace ilang {

/// \brief the class to extract counterexample from a vcd file
class CexExtractor {
public:
  // -------------------- TYPES ------------------ //
  /// a uniform value type (any radix/bool/bv should be fine)
  typedef std::string vlg_val;
  /// val_name -> val
  typedef std::map<std::string, vlg_val>  cex_t;
  /// when generating cex, will only use the regs
  typedef std::map<std::string, bool>  cex_is_reg_t;
  /// a function to determine if some name is a true signal
  /// and a register in the original design or not
  typedef std::function<bool(const std::string &)> is_reg_t;

protected:
  /// the stored cex
  cex_t cex;
  /// whether each var is reg or not
  cex_is_reg_t cex_is_reg;
  /// the helper function to extract info from vcd
  /// for future extension, you can replace this function
  /// to deal with other file format
  void virtual parse_from(const std::string & vcd_file_name, 
    const std::string & scope, is_reg_t is_reg, bool reg_only);

public:
  // -------------------- CONSTRUCTOR ------------------ //
  /// \brief to specify the input vcd name
  /// and also the scope name (the submodule instance name)
  /// to look at
  CexExtractor(const std::string & vcd_file_name, 
    const std::string & scope, is_reg_t is_reg, bool reg_only);

  // -------------------- MEMBERS ------------------ //
  /// return a string to be added to the design
  /// the argument actually has no use at all
  std::string GenInvAssert(const std::string & prefix) const;
  /// allow direct access to the counterexample
  const cex_t & GetCex() const;

}; // class CexExtractor


}; // namespace ilang

#endif // CEX_EXTRACT_H__