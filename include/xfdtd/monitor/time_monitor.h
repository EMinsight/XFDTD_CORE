#ifndef __XFDTD_CORE_TIME_MONITOR_H__
#define __XFDTD_CORE_TIME_MONITOR_H__

#include <xfdtd/monitor/monitor.h>

#include <memory>
#include <string>

namespace xfdtd {

class TimeMonitor : public Monitor {
 public:
  TimeMonitor(std::string name, std::unique_ptr<Cube> cube,
              std::string output_dir);

  auto initTimeDependentVariable() -> void override = 0;

  auto output() -> void override;

  auto toString() const -> std::string override;

  auto time() const -> const Array1D<Real>&;

 protected:
  auto time() -> Array1D<Real>&;

  auto setTime(Array1D<Real> time) -> void;

 private:
  Array1D<Real> _time;
};

}  // namespace xfdtd

#endif  // __XFDTD_CORE_TIME_MONITOR_H__
