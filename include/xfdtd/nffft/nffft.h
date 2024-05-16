#ifndef _XFDTD_CORE_NFFFT_H_
#define _XFDTD_CORE_NFFFT_H_

#include <xfdtd/calculation_param/calculation_param.h>
#include <xfdtd/common/index_task.h>
#include <xfdtd/common/type_define.h>
#include <xfdtd/coordinate_system/coordinate_system.h>
#include <xfdtd/electromagnetic_field/electromagnetic_field.h>
#include <xfdtd/grid_space/grid_space.h>
#include <xfdtd/parallel/mpi_config.h>

#include <memory>
#include <string>

namespace xfdtd {

class XFDTDNFFFTException : public XFDTDException {
 public:
  explicit XFDTDNFFFTException(const std::string& message)
      : XFDTDException(message) {}
};

class NFFFT {
 public:
  NFFFT(Index distance_x, Index distance_y, Index distance_z);

  NFFFT(Index distance_x, Index distance_y, Index distance_z,
        std::string_view output_dir);

  NFFFT(const NFFFT&) = delete;

  NFFFT(NFFFT&&) noexcept = default;

  NFFFT& operator=(const NFFFT&) = delete;

  NFFFT& operator=(NFFFT&&) noexcept = default;

  virtual ~NFFFT();

  virtual auto init(std::shared_ptr<const GridSpace> grid_space,
                    std::shared_ptr<const CalculationParam> calculation_param,
                    std::shared_ptr<const EMF> emf) -> void = 0;

  auto initParallelizedConfig() -> void;

  virtual auto initTimeDependentVariable() -> void = 0;

  virtual auto update() -> void = 0;

  auto outputDir() const -> std::string;

  auto distanceX() const -> Index;

  auto distanceY() const -> Index;

  auto distanceZ() const -> Index;

  auto globalBox() const -> GridBox;

  auto nodeBox() const -> GridBox;

  auto globalTaskSurfaceXN() const -> IndexTask;

  auto globalTaskSurfaceXP() const -> IndexTask;

  auto globalTaskSurfaceYN() const -> IndexTask;

  auto globalTaskSurfaceYP() const -> IndexTask;

  auto globalTaskSurfaceZN() const -> IndexTask;

  auto globalTaskSurfaceZP() const -> IndexTask;

  auto nodeTaskSurfaceXN() const -> IndexTask;

  auto nodeTaskSurfaceXP() const -> IndexTask;

  auto nodeTaskSurfaceYN() const -> IndexTask;

  auto nodeTaskSurfaceYP() const -> IndexTask;

  auto nodeTaskSurfaceZN() const -> IndexTask;

  auto nodeTaskSurfaceZP() const -> IndexTask;

  auto valid() const -> bool;

  auto nffftMPIConfig() const -> const MpiConfig&;

  auto setOutputDir(std::string_view output_dir) -> void;

 protected:
  auto defaultInit(std::shared_ptr<const GridSpace> grid_space,
                   std::shared_ptr<const CalculationParam> calculation_param,
                   std::shared_ptr<const EMF> emf) -> void;

  auto setGlobalBox(GridBox box) -> void;

  auto setNodeBox(GridBox box) -> void;

  auto setGlobalTaskSurfaceXN(IndexTask task) -> void;

  auto setGlobalTaskSurfaceXP(IndexTask task) -> void;

  auto setGlobalTaskSurfaceYN(IndexTask task) -> void;

  auto setGlobalTaskSurfaceYP(IndexTask task) -> void;

  auto setGlobalTaskSurfaceZN(IndexTask task) -> void;

  auto setGlobalTaskSurfaceZP(IndexTask task) -> void;

  auto setNodeTaskSurfaceXN(IndexTask task) -> void;

  auto setNodeTaskSurfaceXP(IndexTask task) -> void;

  auto setNodeTaskSurfaceYN(IndexTask task) -> void;

  auto setNodeTaskSurfaceYP(IndexTask task) -> void;

  auto setNodeTaskSurfaceZN(IndexTask task) -> void;

  auto setNodeTaskSurfaceZP(IndexTask task) -> void;

  auto gridSpace() const -> std::shared_ptr<const GridSpace>;

  auto calculationParam() const -> std::shared_ptr<const CalculationParam>;

  auto emf() const -> std::shared_ptr<const EMF>;

 private:
  Index _distance_x, _distance_y, _distance_z;

  std::string _output_dir;

  std::shared_ptr<const GridSpace> _grid_space;
  std::shared_ptr<const CalculationParam> _calculation_param;
  std::shared_ptr<const EMF> _emf;

  GridBox _global_box, _node_box;

  IndexTask _global_task_surface_xn, _global_task_surface_xp,
      _global_task_surface_yn, _global_task_surface_yp, _global_task_surface_zn,
      _global_task_surface_zp;

  IndexTask _node_task_surface_xn, _node_task_surface_xp, _node_task_surface_yn,
      _node_task_surface_yp, _node_task_surface_zn, _node_task_surface_zp;

  MpiConfig _nffft_mpi_config;

  auto initGlobal() -> void;

  auto initNode() -> void;

  auto makeNodeAxisTask(const Axis::Direction& direction) -> IndexTask;
};

}  // namespace xfdtd

#endif  // _XFDTD_CORE_NFFFT_H_
