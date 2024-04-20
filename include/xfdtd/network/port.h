#ifndef _XFDTD_CORE_PORT_H_
#define _XFDTD_CORE_PORT_H_

#include <xfdtd/monitor/current_monitor.h>
#include <xfdtd/monitor/voltage_monitor.h>

#include <complex>
#include <cstddef>
#include <memory>

#include "xfdtd/common/type_define.h"

namespace xfdtd {

class Port {
 public:
  inline static constexpr int CURRENT_TAG = 31;
  inline static constexpr int VOLTAGE_TAG = 32;

 public:
  Port(std::size_t index, bool is_source, std::complex<Real> impedance,
       std::shared_ptr<CurrentMonitor> current_monitor,
       std::shared_ptr<VoltageMonitor> voltage_monitor);

  Port(const Port& port) = delete;

  Port(Port&& port) noexcept = default;

  Port& operator=(const Port& port) = delete;

  Port& operator=(Port&& port) noexcept = default;

  ~Port() = default;

  std::size_t index() const;

  bool isSource() const;

  std::complex<Real> impedance() const;

  const std::shared_ptr<CurrentMonitor>& currentMonitor() const;

  const std::shared_ptr<VoltageMonitor>& voltageMonitor() const;

  const Array1D<std::complex<Real>>& a() const;

  const Array1D<std::complex<Real>>& b() const;

  void calculateSParameters(const Array1D<Real>& frequencies);

  void init(const std::shared_ptr<const GridSpace>& grid_space,
            const std::shared_ptr<const CalculationParam>& calculation_param,
            const std::shared_ptr<const EMF>& emf);

 private:
  std::size_t _index;
  bool _is_source;
  std::complex<Real> _impedance;
  std::shared_ptr<CurrentMonitor> _current_monitor;
  std::shared_ptr<VoltageMonitor> _voltage_monitor;

  Real _dt;

  Array1D<std::complex<Real>> _a, _b;
};

}  // namespace xfdtd

#endif  // _XFDTD_CORE_PORT_H_
