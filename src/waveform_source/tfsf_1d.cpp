#include <xfdtd/common/constant.h>
#include <xfdtd/common/index_task.h>
#include <xfdtd/common/type_define.h>
#include <xfdtd/grid_space/grid_space.h>
#include <xfdtd/waveform_source/tfsf.h>
#include <xfdtd/waveform_source/tfsf_1d.h>
#include <xfdtd/waveform_source/waveform_source.h>

#include <memory>

#include "waveform_source/tfsf_corrector.h"

namespace xfdtd {

TFSF1D::TFSF1D(Index z, bool forward, std::unique_ptr<Waveform> waveform)
    : TFSF{0,
           0,
           z,
           (forward) ? (0) : (constant::PI),
           0,
           0,
           std::move(waveform)},
      _forward{forward},
      _start{z} {}

void TFSF1D::init(std::shared_ptr<GridSpace> grid_space,
                  std::shared_ptr<CalculationParam> calculation_param,
                  std::shared_ptr<EMF> emf) {
  defaultInit(std::move(grid_space), std::move(calculation_param),
              std::move(emf));
}

auto TFSF1D::generateCorrector(const IndexTask &task)
    -> std::unique_ptr<Corrector> {
  if (!task.valid()) {
    return nullptr;
  }
  auto node_global_start =
      task.zRange().start() + gridSpacePtr()->globalBox().origin().i();
  auto node_global_end =
      task.zRange().end() + gridSpacePtr()->globalBox().origin().i();

  if (_start <= node_global_start) {
    return nullptr;
  }

  if (node_global_end <= _start) {
    return nullptr;
  }

  return std::make_unique<TFSF1DCorrector>(
      Grid{0, 0, _start}, gridSpace(), calculationParam(), emf(),
      waveform()->value(), _forward, _projection_x_int, _projection_y_int,
      _projection_z_int, _projection_x_half, _projection_y_half,
      _projection_z_half, _ex_inc, _ey_inc, _ez_inc, _hx_inc, _hy_inc, _hz_inc,
      cax(), cbx(), cay(), cby(), caz(), cbz());
}

}  // namespace xfdtd
