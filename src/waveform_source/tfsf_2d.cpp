#include <xfdtd/common/constant.h>
#include <xfdtd/common/index_task.h>
#include <xfdtd/waveform_source/tfsf_2d.h>

#include <memory>

#include "waveform_source/tfsf_corrector.h"

namespace xfdtd {

TFSF2D::TFSF2D(std::size_t distance_x, std::size_t distance_y, Real phi,
               std::unique_ptr<Waveform> waveform)
    : TFSF{distance_x, distance_y,         0, constant::PI * 0.5, phi,
           0,          std::move(waveform)} {}

void TFSF2D::init(std::shared_ptr<GridSpace> grid_space,
                  std::shared_ptr<CalculationParam> calculation_param,
                  std::shared_ptr<EMF> emf) {
  defaultInit(std::move(grid_space), std::move(calculation_param),
              std::move(emf));
}

std::unique_ptr<Corrector> TFSF2D::generateCorrector(const IndexTask &task) {
  if (!task.valid()) {
    return nullptr;
  }

  // node task to global task
  auto offset = gridSpace()->globalBox().origin();

  auto node_task =
      makeIndexTask(task.xRange() + offset.i(), task.yRange() + offset.j(),
                    task.zRange() + offset.k());
  auto g_task = globalTask();

  auto intersection_task = taskIntersection(node_task, g_task);

  if (!intersection_task.has_value()) {
    return nullptr;
  }

  return std::make_unique<TFSF2DCorrector>(
      intersection_task.value(), nodeTask(), globalTask(),
      gridSpace().get(), calculationParam().get(), emf().get(),
      &waveform()->value(), gridSpace()->globalBox().origin().i(),
      gridSpace()->globalBox().origin().j(), &_projection_x_int,
      &_projection_y_int, &_projection_z_int, &_projection_x_half,
      &_projection_y_half, &_projection_z_half, &_e_inc, &_h_inc, cax(), cbx(),
      cay(), cby(), caz(), cbz(), _transform_e.x(), _transform_e.y(),
      _transform_e.z(), _transform_h.x(), _transform_h.y(), _transform_h.z());
}

}  // namespace xfdtd
