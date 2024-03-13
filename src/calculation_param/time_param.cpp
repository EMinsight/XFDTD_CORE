#include <xfdtd/calculation_param/calculation_param.h>

#include "xfdtd/util/constant.h"

namespace xfdtd {

TimeParam::TimeParam(double cfl) : _cfl{cfl} {}

double TimeParam::calculateDt(double cfl, double dx, double dy, double dz) {
  return cfl / (constant::C_0 *
                std::sqrt(1.0 / (dx * dx) + 1.0 / (dy * dy) + 1.0 / (dz * dz)));
}

double TimeParam::calculateDt(double cfl, double dx, double dy) {
  return cfl / (constant::C_0 * std::sqrt(1.0 / (dx * dx) + 1.0 / (dy * dy)));
}

double TimeParam::calculateDt(double cfl, double dz) {
  return cfl / (constant::C_0 * std::sqrt(1.0 / (dz * dz)));
}

TimeParam::TimeParam(double dt, std::size_t size, std::size_t start_time_step)
    : _dt{dt},
      _start_time_step{start_time_step},
      _size{size},
      _current_time_step{start_time_step} {}

double TimeParam::dt() const { return _dt; }

double TimeParam::cfl() const { return _cfl; }

std::size_t TimeParam::startTimeStep() const { return _start_time_step; }

std::size_t TimeParam::size() const { return _size; }

std::size_t TimeParam::endTimeStep() const { return _start_time_step + _size; }

std::size_t TimeParam::currentTimeStep() const { return _current_time_step; }

std::size_t TimeParam::remainingTimeStep() const {
  return _start_time_step + _size - _current_time_step;
}

void TimeParam::nextStep() { ++_current_time_step; }

void TimeParam::reset() { _current_time_step = _start_time_step; }

void TimeParam::setDt(double dt) { _dt = dt; }

void TimeParam::setTimeParamRunRange(std::size_t end_time_step,
                                     std::size_t start_time_step) {
  if (end_time_step < start_time_step) {
    throw XFDTDException("end_time_step must be greater than start_time_step");
  }

  _start_time_step = start_time_step;
  _size = end_time_step - start_time_step;
  _current_time_step = start_time_step;
}

xt::xarray<double> TimeParam::eTime() const {
  // return xt::arange<double>(_start_time_step + 1, _start_time_step + _size +
  // 1) * _dt;
  return xt::linspace((_start_time_step + 1) * _dt,
                      (_start_time_step + _size) * _dt, _size);
}

xt::xarray<double> TimeParam::hTime() const {
  // return xt::arange<double>(_start_time_step + 0.5,
  //                           _start_time_step + _size + 0.5) *
  //        _dt;
  return xt::linspace((_start_time_step + 0.5) * _dt,
                      (_start_time_step + _size + 0.5) * _dt, _size);
}

}  // namespace xfdtd
