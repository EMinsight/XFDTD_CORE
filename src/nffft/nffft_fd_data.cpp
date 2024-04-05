#include "nffft/nffft_fd_data.h"

#include <xfdtd/nffft/nffft.h>
#include <xfdtd/util/constant.h>
#include <xfdtd/util/transform.h>

#include "xfdtd/coordinate_system/coordinate_system.h"

namespace xfdtd {

FDPlaneData::FDPlaneData(std::shared_ptr<const GridSpace> grid_space,
                         std::shared_ptr<const EMF> emf, double freq,
                         const Divider::IndexTask& task_xn,
                         const Divider::IndexTask& task_xp,
                         const Divider::IndexTask& task_yn,
                         const Divider::IndexTask& task_yp,
                         const Divider::IndexTask& task_zn,
                         const Divider::IndexTask& task_zp)
    : _grid_space{std::move(grid_space)},
      _emf{std::move(emf)},
      _freq{freq},
      _task_xn{task_xn},
      _task_xp{task_xp},
      _task_yn{task_yn},
      _task_yp{task_yp},
      _task_zn{task_zn},
      _task_zp{task_zp} {
  auto generate_surface = [](const Axis::Direction& direction,
                             const Divider::IndexTask& task, auto& ja, auto& jb,
                             auto& ma, auto& mb) {
    if (!task.valid()) {
      return;
    }

    const auto [range_a, range_b, range_c] = transform::xYZToABC(
        std::tuple{task.xRange(), task.yRange(), task.zRange()},
        Axis::formDirectionToXYZ(direction));

    if (range_c.size() != 1) {
      throw XFDTDNFFFTException("FDPlaneData: Major axis size != 1");
    }

    ja = xt::zeros<std::complex<double>>(
        {task.xRange().size(), task.yRange().size(), task.zRange().size()});
    jb = xt::zeros<std::complex<double>>(
        {task.xRange().size(), task.yRange().size(), task.zRange().size()});
    ma = xt::zeros<std::complex<double>>(
        {task.xRange().size(), task.yRange().size(), task.zRange().size()});
    mb = xt::zeros<std::complex<double>>(
        {task.xRange().size(), task.yRange().size(), task.zRange().size()});
  };

  generate_surface(Axis::Direction::XN, task_xn, _jy_xn, _jz_xn, _my_xn,
                   _mz_xn);
  generate_surface(Axis::Direction::XP, task_xp, _jy_xp, _jz_xp, _my_xp,
                   _mz_xp);
  generate_surface(Axis::Direction::YN, task_yn, _jz_yn, _jx_yn, _mz_yn,
                   _mx_yn);
  generate_surface(Axis::Direction::YP, task_yp, _jz_yp, _jx_yp, _mz_yp,
                   _mx_yp);
  generate_surface(Axis::Direction::ZN, task_zn, _jx_zn, _jy_zn, _mx_zn,
                   _my_zn);
  generate_surface(Axis::Direction::ZP, task_zp, _jx_zp, _jy_zp, _mx_zp,
                   _my_zp);
}

auto FDPlaneData::initDFT(std::size_t total_time_step, double dt) -> void {
  using namespace std::complex_literals;
  _transform_e = xt::zeros<std::complex<double>>({total_time_step});
  _transform_h = xt::zeros<std::complex<double>>({total_time_step});

  for (std::size_t t{0}; t < total_time_step; ++t) {
    _transform_e(t) =
        dt * std::exp(-1i * 2.0 * constant::PI * (_freq * (t + 1) * dt));
    _transform_h(t) =
        dt * std::exp(-1i * 2.0 * constant::PI * (_freq * (t + 0.5) * dt));
  }
}

auto FDPlaneData::frequency() const -> double { return _freq; }

}  // namespace xfdtd
