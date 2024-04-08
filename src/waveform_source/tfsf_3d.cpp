#include <xfdtd/waveform_source/tfsf_3d.h>
#include <xfdtd/common/index_task.h>
#include <utility>

#include "waveform_source/tfsf_corrector.h"

namespace xfdtd {

TFSF3D::TFSF3D(std::size_t x, std::size_t y, std::size_t z, Real theta,
               Real phi, Real psi, std::unique_ptr<Waveform> waveform)
    : TFSF{x, y, z, theta, phi, psi, std::move(waveform)} {}

void TFSF3D::init(std::shared_ptr<GridSpace> grid_space,
                  std::shared_ptr<CalculationParam> calculation_param,
                  std::shared_ptr<EMF> emf) {
  defaultInit(std::move(grid_space), std::move(calculation_param),
              std::move(emf));
}

std::unique_ptr<Corrector> TFSF3D::generateCorrector(
    const IndexTask& task) {
  auto global_ey_task_xn = nodeEyTaskXN(task);

  auto global_ez_task_xn = nodeEzTaskXN(task);

  auto global_ey_task_xp = nodeEyTaskXP(task);

  auto global_ez_task_xp = nodeEzTaskXP(task);

  auto global_ex_task_yn = nodeExTaskYN(task);

  auto global_ez_task_yn = nodeEzTaskYN(task);

  auto global_ex_task_yp = nodeExTaskYP(task);

  auto global_ez_task_yp = nodeEzTaskYP(task);

  auto global_ex_task_zn = nodeExTaskZN(task);

  auto global_ey_task_zn = nodeEyTaskZN(task);

  auto global_ex_task_zp = nodeExTaskZP(task);

  auto global_ey_task_zp = nodeEyTaskZP(task);

  if (!global_ey_task_xn.valid() && !global_ez_task_xn.valid() &&
      !global_ey_task_xp.valid() && !global_ez_task_xp.valid() &&
      !global_ex_task_yn.valid() && !global_ez_task_yn.valid() &&
      !global_ex_task_yp.valid() && !global_ez_task_yp.valid() &&
      !global_ex_task_zn.valid() && !global_ey_task_zn.valid() &&
      !global_ex_task_zp.valid() && !global_ey_task_zp.valid()) {
    return nullptr;
  }

  return std::make_unique<TFSF3DCorrector>(
      globalBox().origin(), gridSpace(), calculationParam(), emf(),
      waveform()->value(), global_ey_task_xn, global_ez_task_xn,
      global_ey_task_xp, global_ez_task_xp, global_ez_task_yn,
      global_ex_task_yn, global_ez_task_yp, global_ex_task_yp,
      global_ex_task_zn, global_ey_task_zn, global_ex_task_zp,
      global_ey_task_zp, _projection_x_int, _projection_y_int,
      _projection_z_int, _projection_x_half, _projection_y_half,
      _projection_z_half, _ex_inc, _ey_inc, _ez_inc, _hx_inc, _hy_inc, _hz_inc,
      cax(), cbx(), cay(), cby(), caz(), cbz());
}

}  // namespace xfdtd
