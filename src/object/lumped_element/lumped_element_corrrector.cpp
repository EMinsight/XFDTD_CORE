#include <xtensor.hpp>
#include <xtensor/xview.hpp>

#include "object/lumped_element_corrector.h"

namespace xfdtd {

void VoltageSourceCorrector::correctE() {
  auto&& coeff_view = xt::view(
      _coeff_v, xt::range(_local_task._x_range[0], _local_task._x_range[1]),
      xt::range(_local_task._y_range[0], _local_task._y_range[1]),
      xt::range(_local_task._z_range[0], _local_task._z_range[1]));

  auto&& e_view =
      xt::view(_e_field, xt::range(_task._x_range[0], _task._x_range[1]),
               xt::range(_task._y_range[0], _task._y_range[1]),
               xt::range(_task._z_range[0], _task._z_range[1]));
  e_view += coeff_view *
            _waveform(_calculation_param->timeParam()->currentTimeStep());
}

void VoltageSourceCorrector::correctH() {}

void CurrentSourceCorrector::correctE() {
  auto&& coeff_view = xt::view(
      _coeff_i, xt::range(_local_task._x_range[0], _local_task._x_range[1]),
      xt::range(_local_task._y_range[0], _local_task._y_range[1]),
      xt::range(_local_task._z_range[0], _local_task._z_range[1]));

  auto&& e_view =
      xt::view(_e_field, xt::range(_task._x_range[0], _task._x_range[1]),
               xt::range(_task._y_range[0], _task._y_range[1]),
               xt::range(_task._z_range[0], _task._z_range[1]));
  e_view += coeff_view *
            _waveform(_calculation_param->timeParam()->currentTimeStep());
}

void CurrentSourceCorrector::correctH() {}

void InductorCorrector::correctE() {
  //   auto&& coeff_ecj = xt::view(
  //       _cecjc, xt::range(_local_task._x_range[0], _local_task._x_range[1]),
  //       xt::range(_local_task._y_range[0], _local_task._y_range[1]),
  //       xt::range(_local_task._z_range[0], _local_task._z_range[1]));
  //   auto&& coeff_jce = xt::view(
  //       _cjcec, xt::range(_local_task._x_range[0], _local_task._x_range[1]),
  //       xt::range(_local_task._y_range[0], _local_task._y_range[1]),
  //       xt::range(_local_task._z_range[0], _local_task._z_range[1]));
  //   auto&& j_view =
  //       xt::view(_j, xt::range(_local_task._x_range[0],
  //       _local_task._x_range[1]),
  //                xt::range(_local_task._y_range[0], _local_task._y_range[1]),
  //                xt::range(_local_task._z_range[0],
  //                _local_task._z_range[1]));

  //   auto&& e_view =
  //       xt::view(_e_field, xt::range(_task._x_range[0], _task._x_range[1]),
  //                xt::range(_task._y_range[0], _task._y_range[1]),
  //                xt::range(_task._z_range[0], _task._z_range[1]));

  //   e_view += coeff_ecj * j_view;
  //   j_view += coeff_jce * e_view;
  const auto local_is = _local_task._x_range[0];
  const auto local_ie = _local_task._x_range[1];
  const auto local_js = _local_task._y_range[0];
  const auto local_je = _local_task._y_range[1];
  const auto local_ks = _local_task._z_range[0];
  const auto local_ke = _local_task._z_range[1];

  const auto is = _task._x_range[0];
  const auto ie = _task._x_range[1];
  const auto js = _task._y_range[0];
  const auto je = _task._y_range[1];
  const auto ks = _task._z_range[0];
  const auto ke = _task._z_range[1];

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
