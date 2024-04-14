#ifndef _XFDTD_CORE_FIELD_MONITOR_H_
#define _XFDTD_CORE_FIELD_MONITOR_H_

#include <xfdtd/electromagnetic_field/electromagnetic_field.h>
#include <xfdtd/monitor/monitor.h>

#include <memory>

namespace xfdtd {

class FieldMonitor : public Monitor {
 public:
  FieldMonitor(std::unique_ptr<Shape> shape, EMF::Field field,
               std::string name = "feild_monitor",
               std::string output_dir_path = "xfdtd_output");

  FieldMonitor(const FieldMonitor&) = delete;

  FieldMonitor(FieldMonitor&&) noexcept = default;

  FieldMonitor& operator=(const FieldMonitor&) = delete;

  FieldMonitor& operator=(FieldMonitor&&) noexcept = default;

  ~FieldMonitor() override = default;

  void init(std::shared_ptr<const GridSpace> grid_space,
            std::shared_ptr<const CalculationParam> calculation_param,
            std::shared_ptr<const EMF> emf) override;

  void update() override;

  void output() override;

  Axis::XYZ axis() const;

  EMF::Field field() const;

  auto initParallelizedConfig() -> void override;

 protected:
  auto gatherData() -> void override;

 private:
  EMF::Field _field;

  std::vector<MpiSupport::Block::Profile> _profiles;
  std::vector<MpiSupport::Block> _blocks_mpi;
};

}  // namespace xfdtd

#endif  // _XFDTD_CORE_FIELD_MONITOR_H_
