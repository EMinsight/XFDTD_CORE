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

  auto g_task = globalTask();

  auto intersection_task = taskIntersection(task, g_task);

  if (!intersection_task.has_value()) {
    return nullptr;
  }

  return std::make_unique<TFSF1DCorrector>(
      intersection_task.value(), gridSpace().get(), calculationParam().get(),
      emf().get(), &waveform()->value(), globalTask(),
      gridSpace()->globalBox().origin().k(), &_projection_x_int,
      &_projection_y_int, &_projection_z_int, &_projection_x_half,
      &_projection_y_half, &_projection_z_half, &_e_inc, &_h_inc, cax(), cbx(),
      cay(), cby(), caz(), cbz(), _transform_e.x(), _transform_e.y(),
      _transform_e.z(), _transform_h.x(), _transform_h.y(), _transform_h.z(),
      _forward);
}

}  // namespace xfdtd
