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

auto TimeMonitor::time() const -> const xt::xarray<double>& { return _time; }

auto TimeMonitor::time() -> xt::xarray<double>& { return _time; }

auto TimeMonitor::setTime(xt::xarray<double> time) -> void {
  _time = std::move(time);
}

}  // namespace xfdtd
