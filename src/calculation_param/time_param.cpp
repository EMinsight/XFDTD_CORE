#include <xfdtd/calculation_param/time_param.h>
#include <xfdtd/common/type_define.h>

#include "xfdtd/common/constant.h"

namespace xfdtd {

TimeParam::TimeParam(Real cfl) : _cfl{cfl} {}

Real TimeParam::calculateDt(Real cfl, Real dx, Real dy, Real dz) {
  return cfl / (constant::C_0 *
                std::sqrt(1.0 / (dx * dx) + 1.0 / (dy * dy) + 1.0 / (dz * dz)));
}

Real TimeParam::calculateDt(Real cfl, Real dx, Real dy) {
  return cfl / (constant::C_0 * std::sqrt(1.0 / (dx * dx) + 1.0 / (dy * dy)));
}

Real TimeParam::calculateDt(Real cfl, Real dz) {
  return cfl / (constant::C_0 * std::sqrt(1.0 / (dz * dz)));
}

TimeParam::TimeParam(Real dt, std::size_t size, std::size_t start_time_step)
    : _dt{dt},
      _start_time_step{start_time_step},
      _size{size},
      _current_time_step{start_time_step} {}

Real TimeParam::dt() const { return _dt; }

Real TimeParam::cfl() const { return _cfl; }

std::size_t TimeParam::startTimeStep() const { return _start_time_step; }

std::size_t TimeParam::size() const { return _size; }

std::size_t TimeParam::endTimeStep() const { return _start_time_step + _size; }

std::size_t TimeParam::currentTimeStep() const { return _current_time_step; }

std::size_t TimeParam::remainingTimeStep() const {
  return _start_time_step + _size - _current_time_step;
}

void TimeParam::nextStep() { ++_current_time_step; }

void TimeParam::reset() { _current_time_step = _start_time_step; }

void TimeParam::setDt(Real dt) { _dt = dt; }

void TimeParam::setTimeParamRunRange(std::size_t end_time_step,
                                     std::size_t start_time_step) {
  if (end_time_step < start_time_step) {
    throw XFDTDTimeParamException(
        "end_time_step must be greater than start_time_step");
  }

  _start_time_step = start_time_step;
  _size = end_time_step - start_time_step;
  _current_time_step = start_time_step;
}

auto TimeParam::eTime() const -> Array1D<Real> {
  return xt::linspace((_start_time_step + 1) * _dt,
                      (_start_time_step + _size) * _dt, _size);
}

auto TimeParam::hTime() const -> Array1D<Real> {
  return xt::linspace((_start_time_step + 0.5) * _dt,
                      (_start_time_step + _size + 0.5) * _dt, _size);
}

}  // namespace xfdtd
