#include <xfdtd/network/port.h>

#include <complex>

#include "xfdtd/util/dft.h"

namespace xfdtd {

Port::Port(std::size_t index, bool is_source, std::complex<double> impedance,
           std::shared_ptr<CurrentMonitor> current_monitor,
           std::shared_ptr<VoltageMonitor> voltage_monitor)
    : _index{index},
      _is_source{is_source},
      _impedance{impedance},
      _current_monitor{std::move(current_monitor)},
      _voltage_monitor{std::move(voltage_monitor)} {}

std::size_t Port::index() const { return _index; }

bool Port::isSource() const { return _is_source; }

std::complex<double> Port::impedance() const { return _impedance; }

const std::shared_ptr<CurrentMonitor>& Port::currentMonitor() const {
  return _current_monitor;
  }

const std::shared_ptr<VoltageMonitor>& Port::voltageMonitor() const {
  return _voltage_monitor;
}

const xt::xarray<std::complex<double>>& Port::a() const { return _a; }

const xt::xarray<std::complex<double>>& Port::b() const { return _b; }

void Port::init(
    const std::shared_ptr<const GridSpace>& grid_space,
    const std::shared_ptr<const CalculationParam>& calculation_param,
    const std::shared_ptr<const EMF>& emf) {
  _dt = calculation_param->timeParam()->dt();
}

void Port::calculateSParameters(const xt::xarray<double>& frequencies) {
  auto current = _current_monitor->data();
  auto voltage = _voltage_monitor->data();

  auto dt{_dt};
  auto z{_impedance};
  auto k{std::sqrt(std::real(z))};
  auto i{dft(current, dt, frequencies, -0.5 * dt)};
  auto v{dft(voltage, dt, frequencies)};

  _a = 0.5 * (v + z * i) / k;
  _b = 0.5 * (v - std::conj(z) * i) / k;
}

}  // namespace xfdtd
