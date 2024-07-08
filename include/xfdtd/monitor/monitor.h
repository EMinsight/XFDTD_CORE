#ifndef _XFDTD_CORE_MONITOR_H_
#define _XFDTD_CORE_MONITOR_H_

#include <xfdtd/calculation_param/calculation_param.h>
#include <xfdtd/common/index_task.h>
#include <xfdtd/electromagnetic_field/electromagnetic_field.h>
#include <xfdtd/exception/exception.h>
#include <xfdtd/grid_space/grid_space.h>
#include <xfdtd/parallel/mpi_config.h>
#include <xfdtd/parallel/mpi_support.h>
#include <xfdtd/parallel/parallelized_config.h>
#include <xfdtd/shape/shape.h>

#include <memory>
#include <string>

namespace xfdtd {

class XFDTDMonitorException : public XFDTDException {
 public:
  explicit XFDTDMonitorException(
      std::string message = "XFDTD Monitor Exception")
      : XFDTDException(std::move(message)) {}
};

class Monitor {
 public:
  explicit Monitor(std::unique_ptr<Shape> shape, std::string name = "monitor",
                   std::string output_dir = "xfdtd_output");

  Monitor(const Monitor&) = delete;

  Monitor(Monitor&&) noexcept = default;

  Monitor& operator=(const Monitor&) = delete;

  Monitor& operator=(Monitor&&) noexcept = default;

  virtual ~Monitor() = default;

  virtual void init(std::shared_ptr<const GridSpace> grid_space,
                    std::shared_ptr<const CalculationParam> calculation_param,
                    std::shared_ptr<const EMF> emf) = 0;

  virtual void update() = 0;

  virtual auto initParallelizedConfig() -> void = 0;

  const std::unique_ptr<Shape>& shape() const;

  const std::string& name() const;

  const std::string& outputDir() const;

  const Array<Real>& data() const;

  std::unique_ptr<Shape>& shape();

  Array<Real>& data();

  void setName(std::string name);

  void setOutputDir(std::string output_dir);

  virtual void output();

  virtual void initTimeDependentVariable();

  GridBox globalGridBox() const;

  GridBox nodeGridBox() const;

  const MpiConfig& monitorMpiConfig() const;

  virtual std::string toString() const;

  auto globalTask() const -> IndexTask;

  auto nodeTask() const -> IndexTask;

  virtual auto valid() const -> bool;

  auto gridSpace() const -> std::shared_ptr<const GridSpace> {
    return _grid_space;
  }

  auto calculationParam() const -> std::shared_ptr<const CalculationParam> {
    return _calculation_param;
  }

  auto emf() const -> std::shared_ptr<const EMF> { return _emf; }

 protected:
  void defaultInit(std::shared_ptr<const GridSpace> grid_space,
                   std::shared_ptr<const CalculationParam> calculation_param,
                   std::shared_ptr<const EMF> emf);

  const GridSpace* gridSpacePtr() const;

  const CalculationParam* calculationParamPtr() const;

  const EMF* emfPtr() const;

  auto nodeGridBox() -> GridBox&;

  auto nodeTask() -> IndexTask&;

  MpiConfig& monitorMpiConfig();

  virtual auto gatherData() -> void;

  auto setGlobalGridBox(GridBox grid_box) -> void;

  auto setNodeGridBox(GridBox grid_box) -> void;

  auto setGlobalTask(IndexTask task) -> void;

  auto setNodeTask(IndexTask task) -> void;

  auto makeMpiSubComm() -> void;

 private:
  std::unique_ptr<Shape> _shape;
  std::string _name;
  std::string _output_dir;
  Array<Real> _data;

  std::shared_ptr<const GridSpace> _grid_space;
  std::shared_ptr<const CalculationParam> _calculation_param;
  std::shared_ptr<const EMF> _emf;

  GridBox _global_grid_box;
  GridBox _node_grid_box;
  IndexTask _global_task;
  IndexTask _node_task;

  MpiConfig _monitor_mpi_config;
};

}  // namespace xfdtd

#endif  // _XFDTD_CORE_MONITOR_H_
