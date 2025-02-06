#ifndef __XFDTD_CORE_SIMULATION_FLAG_H__
#define __XFDTD_CORE_SIMULATION_FLAG_H__

#include <xfdtd/common/type_define.h>

namespace xfdtd {

enum class SimulationInitFlag {
  SimulationStart,
  SimulationEnd,
  InitStart,
  InitEnd,
  UpdateStart,
  UpdateEnd
};

enum class SimulationIteratorFlag {
  UpdateHSStart,
  UpdateHSEnd,
  UpdateEStart,
  UpdateEEnd,
  NextStep
};

class SimulationFlagVisitor {
 public:
  virtual auto initStep(SimulationInitFlag flag) -> void = 0;

  virtual auto iteratorStep(SimulationIteratorFlag flag, Index cur, Index start,
                            Index end) -> void = 0;
};

}  // namespace xfdtd

#endif  // __XFDTD_CORE_SIMULATION_FLAG_H__
