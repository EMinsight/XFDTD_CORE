#include <xtensor.hpp>
#include <xtensor/xview.hpp>

#include "object/lumped_element_corrector.h"

namespace xfdtd {

void VoltageSourceCorrector::correctE() {
  auto&& coeff_view = xt::view(
      _coeff_v, xt::range(_local_task.xRange().start(), _local_task.xRange().end()),
      xt::range(_local_task.yRange().start(), _local_task.yRange().end()),
      xt::range(_local_task.zRange().start(), _local_task.zRange().end()));

  auto&& e_view =
      xt::view(_e_field, xt::range(_task.xRange().start(), _task.xRange().end()),
               xt::range(_task.yRange().start(), _task.yRange().end()),
               xt::range(_task.zRange().start(), _task.zRange().end()));
  e_view += coeff_view *
            _waveform(_calculation_param->timeParam()->currentTimeStep());
}

void VoltageSourceCorrector::correctH() {}

void CurrentSourceCorrector::correctE() {
  auto&& coeff_view = xt::view(
      _coeff_i, xt::range(_local_task.xRange().start(), _local_task.xRange().end()),
      xt::range(_local_task.yRange().start(), _local_task.yRange().end()),
      xt::range(_local_task.zRange().start(), _local_task.zRange().end()));

  auto&& e_view =
      xt::view(_e_field, xt::range(_task.xRange().start(), _task.xRange().end()),
               xt::range(_task.yRange().start(), _task.yRange().end()),
               xt::range(_task.zRange().start(), _task.zRange().end()));
  e_view += coeff_view *
            _waveform(_calculation_param->timeParam()->currentTimeStep());
}

void CurrentSourceCorrector::correctH() {}

void InductorCorrector::correctE() {
  //   auto&& coeff_ecj = xt::view(
  //       _cecjc, xt::range(_local_task.xRange().start(), _local_task.xRange().end()),
  //       xt::range(_local_task.yRange().start(), _local_task.yRange().end()),
  //       xt::range(_local_task.zRange().start(), _local_task.zRange().end()));
  //   auto&& coeff_jce = xt::view(
  //       _cjcec, xt::range(_local_task.xRange().start(), _local_task.xRange().end()),
  //       xt::range(_local_task.yRange().start(), _local_task.yRange().end()),
  //       xt::range(_local_task.zRange().start(), _local_task.zRange().end()));
  //   auto&& j_view =
  //       xt::view(_j, xt::range(_local_task.xRange().start(),
  //       _local_task.xRange().end()),
  //                xt::range(_local_task.yRange().start(), _local_task.yRange().end()),
  //                xt::range(_local_task.zRange().start(),
  //                _local_task.zRange().end()));

  //   auto&& e_view =
  //       xt::view(_e_field, xt::range(_task.xRange().start(), _task.xRange().end()),
  //                xt::range(_task.yRange().start(), _task.yRange().end()),
  //                xt::range(_task.zRange().start(), _task.zRange().end()));

  //   e_view += coeff_ecj * j_view;
  //   j_view += coeff_jce * e_view;
  const auto local_is = _local_task.xRange().start();
  const auto local_ie = _local_task.xRange().end();
  const auto local_js = _local_task.yRange().start();
  const auto local_je = _local_task.yRange().end();
  const auto local_ks = _local_task.zRange().start();
  const auto local_ke = _local_task.zRange().end();

  const auto is = _task.xRange().start();
  const auto ie = _task.xRange().end();
  const auto js = _task.yRange().start();
  const auto je = _task.yRange().end();
  const auto ks = _task.zRange().start();
  const auto ke = _task.zRange().end();

  const auto nx = ie - is;
  const auto ny = je - js;
  const auto nz = ke - ks;

  for (std::size_t i = 0; i < nx; ++i) {
    for (std::size_t j = 0; j < ny; ++j) {
      for (std::size_t k = 0; k < nz; ++k) {
        _e_field(i + is, j + js, k + ks) +=
            _cecjc(i + local_is, j + local_js, k + local_ks) *
            _j(i + local_is, j + local_js, k + local_ks);
        _j(i + local_is, j + local_js, k + local_ks) +=
            _cjcec(i + local_is, j + local_js, k + local_ks) *
            _e_field(i + is, j + js, k + ks);
      }
    }
  }
}

void InductorCorrector::correctH() {}

}  // namespace xfdtd
