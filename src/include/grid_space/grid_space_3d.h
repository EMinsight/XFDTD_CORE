#ifndef _XFDTD_CORE_GRID_SPACE_3D_H_
#define _XFDTD_CORE_GRID_SPACE_3D_H_

#include <xfdtd/grid_space/grid_space.h>

namespace xfdtd {

class GridSpace3D : public GridSpace {
 public:
  GridSpace3D(Real based_dx, Real based_dy, Real based_dz,
              Array1D<Real> e_node_x, Array1D<Real> e_node_y,
              Array1D<Real> e_node_z);

  GridSpace3D(Type type, GridBox global_box, Real based_dx, Real based_dy,
              Real based_dz, Real min_dx, Real min_dy, Real min_dz,
              Array1D<Real> e_node_x, Array1D<Real> e_node_y,
              Array1D<Real> e_node_z, Array1D<Real> h_node_x,
              Array1D<Real> h_node_y, Array1D<Real> h_node_z,
              Array1D<Real> e_size_x, Array1D<Real> e_size_y,
              Array1D<Real> e_size_z, Array1D<Real> h_size_x,
              Array1D<Real> h_size_y, Array1D<Real> h_size_z);

  std::unique_ptr<GridSpace> subGridSpace(std::size_t start_i,
                                          std::size_t start_j,
                                          std::size_t start_k,
                                          std::size_t end_i, std::size_t end_j,
                                          std::size_t end_k) const override;

  void correctGridSpace() override;
};

}  // namespace xfdtd

#endif  // _XFDTD_CORE_GRID_SPACE_3D_H_
