/// file
/// The implementation for LMAC ILA.

#include "lmac.h"

namespace ilang {

Ila GetLmacIla(const std::string& name) {
  auto lmac = Ila("LMAC");
  DefineArchState(lmac);
  DefineInstruction(lmac);
  return lmac;
}

void DefineArchState(Ila& m) {
  // --- implementation dependent
  auto clk = m.NewBoolInput("clk");         // positive edge
  auto reg_clk = m.NewBoolInput("reg_clk"); // positive edge
  auto rst = m.NewBoolInput("rst");         // low true
  auto phy_rst = m.NewBoolInput("phy_rst"); // low true

  // --- static signals
  auto link_up = m.NewBoolState("link_up");     // 1 indicating ethernet is up
  auto mac_ctrl = m.NewBvInput("mac_ctrl", 32); // configuration control signals
  auto mac_ctrl_1 = m.NewBvInput("mac_ctrl_1", 32); // the second set of config.
  auto mac_rxv_en = m.NewBoolInput("mac_rxv_en");   // 1 for receive path enable

  // --- register interface
  auto host_addr = m.NewBvInput("host_addr", 16);
  // byte address for selecting MMIO registers
  // valid 1 clock cycle before the "reg_read_start" signal

  auto reg_read_start = m.NewBoolInput("reg_read_start");
  // pulse (up for one and exactly one cycle)
  // 1 for starting the read for a register in the LMAC core
  // 0 for idle

  auto read_reg_done_out = m.NewBoolState("read_reg_done_out");
  // indicating the register read is done and its data is available on the
  // MAC_REGDOUT bus
  // pulse (up for one and exactly one cycle, even it is and output)
  // 1 for data is ready for the read to the register
  // 0 for data not available
  // * It is the 5-clock delay of the "reg_read_start" signal

  auto mac_reg_d_out = m.NewBvState("mac_reg_d_out", 32);
  // the 32-bit read data

  // --- MMIO registers
#define MMIO_REG_SIZE 32
  auto fmac_tx_pkt_cnt = m.NewBvState("FMAC_TX_PKT_CNT", MMIO_REG_SIZE);
  auto fmac_rx_pkt_cnt_lo = m.NewBvState("FMAC_RX_PKT_CNT_LO", MMIO_REG_SIZE);
  auto fmac_rx_pkt_cnt_hi = m.NewBvState("FMAC_RX_PKT_CNT_HI", MMIO_REG_SIZE);
  auto fmac_tx_byte_cnt = m.NewBvState("FMAC_TX_BYTE_CNT", MMIO_REG_SIZE);
  auto fmac_rx_byte_cnt_lo = m.NewBvState("FMAC_RX_BYTE_CNT_LO", MMIO_REG_SIZE);
  auto fmac_rx_byte_cnt_hi = m.NewBvState("FMAC_RX_BYTE_CNT_HI", MMIO_REG_SIZE);
  // TODO

  // helpers
  auto read_reg_cycle_cnt = m.NewBvState("lamc_read_reg_cycle_cnt", 8);
  auto read_reg_cache_val =
      m.NewBvState("lmac_read_reg_cache_val", MMIO_REG_SIZE);

  return;
}

// XXX invariant:
// 1. read_reg_done_out & reg_read_start = 0
// 2. reg_read_start should be 1 for only one cycle
// 3. read_reg_done_out -> (lmac_read_reg_cycle_cnt == 0)
//
// XXX does the verification target generation handle edge dependent
// properties? (one and exactly one cycle)

void DefineInstruction(Ila& m) {
  // child ILA for 5 cycle delay
  auto reg_read_child = m.NewChild("reg_read_child");
  {
    auto& c = reg_read_child;

    // valid
    auto valid = c.state("lmac_read_reg_cycle_cnt") != 0;
    reg_read_child.SetValid(valid);

    // fetch
    auto fetch = Concat(c.state("lmac_read_reg_cycle_cnt"),
                        c.state("lmac_read_reg_cache_val"));
    reg_read_child.SetFetch(fetch);

    // instructions for 5 cycle delay
    auto cnt = c.state("lmac_read_reg_cycle_cnt");
    auto out = c.state("read_reg_done_out");

    // cnt: 1 -> 2 -> 3 -> 4 -> 5 -> 0
    // out: 0 -> 0 -> 0 -> 0 -> 1 -> 0

    auto instr_delay_ini = c.NewInstr("reg_read_delay_ini");
    {
      instr_delay_ini.SetDecode((cnt > 0) & (cnt < 4));
      instr_delay_ini.SetUpdate(cnt, cnt + 1);
    }

    auto instr_delay_out = c.NewInstr("reg_read_delay_out");
    {
      instr_delay_out.SetDecode(cnt == 4);
      instr_delay_out.SetUpdate(cnt, cnt + 1);
      instr_delay_out.SetUpdate(out, BoolConst(true));
    }

    auto instr_delay_end = c.NewInstr("reg_read_delay_end");
    {
      instr_delay_end.SetDecode(cnt == 5);
      instr_delay_end.SetUpdate(cnt, BvConst(0, 8));
    }
  }

  auto instr_read_tx_pkt_cnt = m.NewInstr("READ_FMAC_TX_PKT_CNT");
  {
    // decode
    auto decode = (m.input("reg_read_start") == 1) &
                  (m.state("read_reg_done_out") == 0) &
                  (m.input("host_addr") == LMAC_MMIO_OFFSET_FMAC_TX_PKT_CNT);
    instr_read_tx_pkt_cnt.SetDecode(decode);

    // updates
    instr_read_tx_pkt_cnt.SetUpdate(m.state("lmac_read_reg_cache_val"),
                                    m.state("fmac_tx_pkt_cnt"));
    instr_read_tx_pkt_cnt.SetUpdate(m.state("lmac_read_reg_cycle_cnt"),
                                    BvConst(1, 8));
  }

  return;
}

}; // namespace ilang

