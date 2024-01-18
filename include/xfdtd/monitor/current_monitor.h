#ifndef _XFDTD_LIB_CURRENT_MONITOR_H_
#define _XFDTD_LIB_CURRENT_MONITOR_H_

#include "xfdtd/monitor/monitor.h"

namespace xfdtd {

class CurrentMonitor : public Monitor {
 public:
  CurrentMonitor(std::string name, std::unique_ptr<Cube> cube,
                 Axis::Direction direction, std::string output_dir);

  CurrentMonitor(const CurrentMonitor&) = delete;

  CurrentMonitor(CurrentMonitor&&) noexcept = default;

  CurrentMonitor& operator=(const CurrentMonitor&) = delete;

  CurrentMonitor& operator=(CurrentMonitor&&) noexcept = default;

  ~CurrentMonitor() override = default;

  void init(std::shared_ptr<const GridSpace> grid_space,
            std::shared_ptr<const CalculationParam> calculation_param,
            std::shared_ptr<const EMF> emf) override;

  void update() override;

 private:
  Axis::Direction _direction;
  std::size_t _is, _ie, _js, _je, _ks, _ke;
  xt::xarray<double> _da, _db;
  double _positive;
};

}  // namespace xfdtd

#endif  // _XFDTD_LIB_CURRENT_MONITOR_H_
