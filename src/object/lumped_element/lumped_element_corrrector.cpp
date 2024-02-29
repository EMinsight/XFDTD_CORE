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
  e_view +=
      _coeff_v * _waveform(_calculation_param->timeParam()->currentTimeStep());
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
  e_view +=
      _coeff_i * _waveform(_calculation_param->timeParam()->currentTimeStep());
}

void CurrentSourceCorrector::correctH() {}

void InductorCorrector::correctE() {
  auto&& coeff_ecj = xt::view(
      _cecjc, xt::range(_local_task._x_range[0], _local_task._x_range[1]),
      xt::range(_local_task._y_range[0], _local_task._y_range[1]),
      xt::range(_local_task._z_range[0], _local_task._z_range[1]));
  auto&& coeff_jce = xt::view(
      _cjcec, xt::range(_local_task._x_range[0], _local_task._x_range[1]),
      xt::range(_local_task._y_range[0], _local_task._y_range[1]),
      xt::range(_local_task._z_range[0], _local_task._z_range[1]));
  auto&& j_view =
      xt::view(_j, xt::range(_local_task._x_range[0], _local_task._x_range[1]),
               xt::range(_local_task._y_range[0], _local_task._y_range[1]),
               xt::range(_local_task._z_range[0], _local_task._z_range[1]));

  auto&& e_view =
      xt::view(_e_field, xt::range(_task._x_range[0], _task._x_range[1]),
               xt::range(_task._y_range[0], _task._y_range[1]),
               xt::range(_task._z_range[0], _task._z_range[1]));

  e_view += coeff_ecj * j_view;
  j_view += coeff_jce * e_view;
}

void InductorCorrector::correctH() {}

}  // namespace xfdtd
