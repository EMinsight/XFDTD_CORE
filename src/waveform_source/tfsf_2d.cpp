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
  auto xn_task = nodeEzTaskXN(task);
  auto xp_task = nodeEzTaskXP(task);
  auto yn_task = nodeEzTaskYN(task);
  auto yp_task = nodeEzTaskYP(task);

  if (!xn_task.valid() && !xp_task.valid() && !yn_task.valid() &&
      !yp_task.valid()) {
    return nullptr;
  }

  return std::make_unique<TFSF2DCorrector>(
      globalBox().origin(), gridSpace(), calculationParam(), emf(),
      waveform()->value(), xn_task, xp_task, yn_task, yp_task,
      _projection_x_int, _projection_y_int, _projection_z_int,
      _projection_x_half, _projection_y_half, _projection_z_half, _ex_inc,
      _ey_inc, _ez_inc, _hx_inc, _hy_inc, _hz_inc, cax(), cbx(), cay(), cby(),
      caz(), cbz());
}

}  // namespace xfdtd
