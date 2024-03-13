#ifndef _XFDTD_CORE_PORT_H_
#define _XFDTD_CORE_PORT_H_

#include <xfdtd/monitor/current_monitor.h>
#include <xfdtd/monitor/voltage_monitor.h>

#include <complex>
#include <cstddef>
#include <memory>

namespace xfdtd {

class Port {
 public:
  Port(std::size_t index, bool is_source, std::complex<double> impedance,
       std::shared_ptr<CurrentMonitor> current_monitor,
       std::shared_ptr<VoltageMonitor> voltage_monitor);

  Port(const Port& port) = delete;

  Port(Port&& port) noexcept = default;

  Port& operator=(const Port& port) = delete;

  Port& operator=(Port&& port) noexcept = default;

  ~Port() = default;

  std::size_t index() const;

  bool isSource() const;

  std::complex<double> impedance() const;

  const std::shared_ptr<CurrentMonitor>& currentMonitor() const;

  const std::shared_ptr<VoltageMonitor>& voltageMonitor() const;

  const xt::xarray<std::complex<double>>& a() const;

  const xt::xarray<std::complex<double>>& b() const;

  void calculateSParameters(const xt::xarray<double>& frequencies);

  void init(const std::shared_ptr<const GridSpace>& grid_space,
            const std::shared_ptr<const CalculationParam>& calculation_param,
            const std::shared_ptr<const EMF>& emf);

 private:
  std::size_t _index;
  bool _is_source;
  std::complex<double> _impedance;
  std::shared_ptr<CurrentMonitor> _current_monitor;
  std::shared_ptr<VoltageMonitor> _voltage_monitor;

  double _dt;

  xt::xarray<std::complex<double>> _a, _b;
};

}  // namespace xfdtd

#endif  // _XFDTD_CORE_PORT_H_
