#include <xfdtd/monitor/time_monitor.h>

#include <utility>

namespace xfdtd {

TimeMonitor::TimeMonitor(std::string name, std::unique_ptr<Cube> cube,
                         std::string output_dir)
    : Monitor{std::move(cube), std::move(name), std::move(output_dir)} {}

auto TimeMonitor::output() -> void {
  if (!valid()) {
    return;
  }

  gatherData();
  auto temp = data();
  data() = xt::stack(xt::xtuple(time(), data()));
  Monitor::output();
  data() = temp;
}

auto TimeMonitor::toString() const -> std::string { return "TimeMonitor"; }

auto TimeMonitor::time() const -> const Array1D<Real>& { return _time; }

auto TimeMonitor::time() -> Array1D<Real>& { return _time; }

auto TimeMonitor::setTime(Array1D<Real> time) -> void {
  _time = std::move(time);
}

}  // namespace xfdtd
