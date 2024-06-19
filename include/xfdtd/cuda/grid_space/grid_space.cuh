#ifndef __XFDTD_CORE_CUDA_GRID_SPACE_GRID_SPACE_CUH__
#define __XFDTD_CORE_CUDA_GRID_SPACE_GRID_SPACE_CUH__

#include <xfdtd/common/type_define.h>

#include <xfdtd/cuda/common.cuh>
#include <xfdtd/cuda/tensor.cuh>
#include <xfdtd/cuda/tensor_hd.cuh>

#include "xfdtd/grid_space/grid_space.h"

namespace xfdtd {

namespace cuda {

/**
 * @brief Holds the grid space data
 *
 */
class GridSpaceData {
  friend class GridSpaceHD;

 public:
  /**
   * @brief e node means that the corner of the cell
   *
   * @return const Array1D<Real>&
   */
  auto eNodeX() const -> const Array1D<Real>&;

  /**
   * @brief h node means that the center of the cell
   *
   * @return const Array1D<Real>&
   */
  auto hNodeX() const -> const Array1D<Real>&;

  auto eNodeY() const -> const Array1D<Real>&;

  auto eNodeZ() const -> const Array1D<Real>&;

  auto hNodeY() const -> const Array1D<Real>&;

  auto hNodeZ() const -> const Array1D<Real>&;

  // auto eNodeX() -> Array1D<Real>&;

  // auto eNodeY() -> Array1D<Real>&;

  // auto eNodeZ() -> Array1D<Real>&;

  // auto hNodeX() -> Array1D<Real>&;

  // auto hNodeY() -> Array1D<Real>&;

  // auto hNodeZ() -> Array1D<Real>&;

  auto eSizeX() const -> const Array1D<Real>&;

  auto eSizeY() const -> const Array1D<Real>&;

  auto eSizeZ() const -> const Array1D<Real>&;

  auto hSizeX() const -> const Array1D<Real>&;

  auto hSizeY() const -> const Array1D<Real>&;

  auto hSizeZ() const -> const Array1D<Real>&;

  /**
   * @brief The number of cells in x direction
   *
   * @return Index
   */
  auto sizeX() const -> Index;

  /**
   * @brief The number of cells in y direction
   *
   * @return Index
   */
  auto sizeY() const -> Index;

  /**
   * @brief The number of cells in z direction
   *
   * @return Index
   */
  auto sizeZ() const -> Index;

 private:
  Real _based_dx, _based_dy, _based_dz;
  Real _min_dx, _min_dy, _min_dz;

  Array1D<Real>*_e_node_x{}, *_e_node_y, *_e_node_z;
  Array1D<Real>*_h_node_x, *_h_node_y, *_h_node_z;
  Array1D<Real>*_e_size_x, *_e_size_y, *_e_size_z;
  Array1D<Real>*_h_size_x, *_h_size_y, *_h_size_z;
};

class GridSpaceHD {
  using Host = xfdtd::GridSpace;
  using Device = GridSpaceData;

 public:
  GridSpaceHD(const Host* host);

  GridSpaceHD(const GridSpaceHD&) = delete;

  GridSpaceHD(GridSpaceHD&&) = delete;

  ~GridSpaceHD();

  auto operator=(const GridSpaceHD&) -> GridSpaceHD& = delete;

  auto operator=(GridSpaceHD&&) -> GridSpaceHD& = delete;

  auto copyHostToDevice() -> void;

  auto copyDeviceToHost() -> void;

  auto releaseDevice() -> void;

  auto host() { return _host; }

  auto device() { return _device; }

  auto host() const { return _host; }

  auto device() const { return _device; }

 private:
  const Host* _host{};
  Device* _device{};

  TensorHD<Real, 1> _e_node_x_hd, _e_node_y_hd, _e_node_z_hd;
  TensorHD<Real, 1> _h_node_x_hd, _h_node_y_hd, _h_node_z_hd;
  TensorHD<Real, 1> _e_size_x_hd, _e_size_y_hd, _e_size_z_hd;
  TensorHD<Real, 1> _h_size_x_hd, _h_size_y_hd, _h_size_z_hd;
};

}  // namespace cuda

}  // namespace xfdtd

#endif  // __XFDTD_CORE_CUDA_GRID_SPACE_GRID_SPACE_CUH__
