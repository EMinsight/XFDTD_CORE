#ifndef __XFDTD_CORE_FIELD_TIME_MONITOR_H__
#define __XFDTD_CORE_FIELD_TIME_MONITOR_H__

#include <xfdtd/monitor/time_monitor.h>

namespace xfdtd {

class FieldTimeMonitor : public TimeMonitor {
 public:
  FieldTimeMonitor(std::string name, std::unique_ptr<Cube> cube,
                   std::string output_dir, EMF::Field field);

  FieldTimeMonitor(const FieldTimeMonitor&) = delete;

  FieldTimeMonitor(FieldTimeMonitor&&) noexcept = default;

  FieldTimeMonitor& operator=(const FieldTimeMonitor&) = delete;

  FieldTimeMonitor& operator=(FieldTimeMonitor&&) noexcept = default;

  ~FieldTimeMonitor() override = default;

  auto init(std::shared_ptr<const GridSpace> grid_space,
            std::shared_ptr<const CalculationParam> calculation_param,
            std::shared_ptr<const EMF> emf) -> void override;

  auto update() -> void override;

  auto initTimeDependentVariable() -> void override;

  auto initParallelizedConfig() -> void override;

  auto toString() const -> std::string override;

  auto field() const -> EMF::Field;

 private:
  EMF::Field _field;

  void* _no_meaningful_data{nullptr};
};

}  // namespace xfdtd

#endif  // __XFDTD_CORE_FIELD_TIME_MONITOR_H__
