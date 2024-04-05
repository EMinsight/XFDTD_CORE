#ifndef _XFDTD_CORE_VOLTAGE_MONITOR_H_
#define _XFDTD_CORE_VOLTAGE_MONITOR_H_

#include <xfdtd/monitor/monitor.h>
#include <xfdtd/monitor/time_monitor.h>

namespace xfdtd {

class VoltageMonitor : public TimeMonitor {
 public:
  VoltageMonitor(std::string name, std::unique_ptr<Cube> cube,
                 Axis::Direction direction, std::string output_dir);

  VoltageMonitor(const VoltageMonitor&) = delete;

  VoltageMonitor(VoltageMonitor&&) noexcept = default;

  VoltageMonitor& operator=(const VoltageMonitor&) = delete;

  VoltageMonitor& operator=(VoltageMonitor&&) noexcept = default;

  ~VoltageMonitor() override = default;

  void init(std::shared_ptr<const GridSpace> grid_space,
            std::shared_ptr<const CalculationParam> calculation_param,
            std::shared_ptr<const EMF> emf) override;

  void update() override;

  void initTimeDependentVariable() override;

  auto initParallelizedConfig() -> void override;

  auto gatherData() -> void override;

  auto toString() const -> std::string override;

 private:
  Axis::Direction _direction;
  std::size_t _is, _ie, _js, _je, _ks, _ke;
  xt::xarray<double> _dc;
  xt::xarray<double> _coff;

  // xt::xarray<double> _time;
  xt::xarray<double> _node_data;
};

}  // namespace xfdtd

#endif  // _XFDTD_CORE_VOLTAGE_MONITOR_H_
