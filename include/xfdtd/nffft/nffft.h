#ifndef _XFDTD_CORE_NFFFT_H_
#define _XFDTD_CORE_NFFFT_H_

#include <xfdtd/calculation_param/calculation_param.h>
#include <xfdtd/common/index_task.h>
#include <xfdtd/common/type_define.h>
#include <xfdtd/coordinate_system/coordinate_system.h>
#include <xfdtd/electromagnetic_field/electromagnetic_field.h>
#include <xfdtd/grid_space/grid_space.h>
#include <xfdtd/parallel/mpi_config.h>

#include <complex>
#include <cstddef>
#include <memory>
#include <string>
#include <vector>
#include <xtensor/xtensor.hpp>

namespace xfdtd {

class XFDTDNFFFTException : public XFDTDException {
 public:
  explicit XFDTDNFFFTException(const std::string& message)
      : XFDTDException(message) {}
};

class FDPlaneData;
class NFFFT {
 public:
  NFFFT(std::size_t distance_x, std::size_t distance_y, std::size_t distance_z,
        Array1D<Real> frequencies, std::string output_dir);

  NFFFT(const NFFFT&) = delete;

  NFFFT(NFFFT&&) noexcept = default;

  NFFFT& operator=(const NFFFT&) = delete;

  NFFFT& operator=(NFFFT&&) noexcept = default;

  ~NFFFT();

  void init(std::shared_ptr<const GridSpace> grid_space,
            std::shared_ptr<const CalculationParam> calculation_param,
            std::shared_ptr<const EMF> emf);

  auto initParallelizedConfig() -> void;

  void initTimeDependentVariable();

  void update();

  auto processFarField(const xt::xtensor<Real, 1>& theta, Real phi,
                       const std::string& sub_dir,
                       const Vector& origin = Vector{0.0, 0.0, 0.0}) const
      -> void;

  auto processFarField(Real theta, const xt::xtensor<Real, 1>& phi,
                       const std::string& sub_dir,
                       const Vector& origin = Vector{0.0, 0.0, 0.0}) const
      -> void;

  void outputRadiationPower();

  void setOutputDir(const std::string& out_dir);

  auto distanceX() const -> std::size_t;

  auto distanceY() const -> std::size_t;

  auto distanceZ() const -> std::size_t;

  auto cube() const -> const Cube&;

  auto globalBox() const -> const GridBox&;

  auto nodeBox() const -> const GridBox&;

  auto globalTaskSurfaceXN() const -> const IndexTask&;

  auto globalTaskSurfaceXP() const -> const IndexTask&;

  auto globalTaskSurfaceYN() const -> const IndexTask&;

  auto globalTaskSurfaceYP() const -> const IndexTask&;

  auto globalTaskSurfaceZN() const -> const IndexTask&;

  auto globalTaskSurfaceZP() const -> const IndexTask&;

  auto nodeTaskSurfaceXN() const -> const IndexTask&;

  auto nodeTaskSurfaceXP() const -> const IndexTask&;

  auto nodeTaskSurfaceYN() const -> const IndexTask&;

  auto nodeTaskSurfaceYP() const -> const IndexTask&;

  auto nodeTaskSurfaceZN() const -> const IndexTask&;

  auto nodeTaskSurfaceZP() const -> const IndexTask&;

  auto valid() const -> bool;

  auto nffftMPIConfig() const -> const MpiConfig&;

 protected:
  auto processFarField(const xt::xtensor<Real, 1>& theta,
                       const xt::xtensor<Real, 1>& phi,
                       const std::string& sub_dir, const Vector& origin) const
      -> void;

  const GridSpace* gridSpacePtr() const;

  const CalculationParam* calculationParamPtr() const;

  const EMF* emfPtr() const;

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

  auto nffftMPIConfig() -> MpiConfig&;

 private:
  std::size_t _distance_x, _distance_y, _distance_z;
  Array1D<Real> _frequencies;
  std::string _output_dir;

  std::shared_ptr<const GridSpace> _grid_space;
  std::shared_ptr<const CalculationParam> _calculation_param;
  std::shared_ptr<const EMF> _emf;

  std::unique_ptr<Cube> _cube;
  GridBox _global_box, _node_box;

  IndexTask _global_task_surface_xn, _global_task_surface_xp,
      _global_task_surface_yn, _global_task_surface_yp, _global_task_surface_zn,
      _global_task_surface_zp;

  IndexTask _node_task_surface_xn, _node_task_surface_xp, _node_task_surface_yn,
      _node_task_surface_yp, _node_task_surface_zn, _node_task_surface_zp;

  MpiConfig _nffft_mpi_config;

  std::vector<FDPlaneData> _fd_plane_data;

  xt::xarray<std::complex<Real>> _transform_e, _transform_h;

  xt::xarray<std::complex<Real>> _a_theta, _a_phi, _f_theta, _f_phi;

  auto initGlobal() -> void;

  auto initNode() -> void;

  auto makeNodeAxisTask(const Axis::Direction& direction) -> IndexTask;

  auto generateSurface() -> void;
};

}  // namespace xfdtd

#endif  // _XFDTD_CORE_NFFFT_H_
