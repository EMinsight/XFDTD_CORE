#ifndef _XFDTD_CORE_GRID_SPACE_2D_H_
#define _XFDTD_CORE_GRID_SPACE_2D_H_

#include <xfdtd/grid_space/grid_space.h>

namespace xfdtd {

class GridSpace2D : public GridSpace {
 public:
  GridSpace2D(Real based_dx, Real based_dy, Array1D<Real> e_node_x,
              Array1D<Real> e_node_y);

  GridSpace2D(Type type, GridBox global_box, Real based_dx, Real based_dy,
              Real min_dx, Real min_dy, Array1D<Real> e_node_x,
              Array1D<Real> e_node_y, Array1D<Real> h_node_x,
              Array1D<Real> h_node_y, Array1D<Real> e_size_x,
              Array1D<Real> e_size_y, Array1D<Real> h_size_x,
              Array1D<Real> h_size_y);

  void correctGridSpace() override;

  std::unique_ptr<GridSpace> subGridSpace(std::size_t start_i,
                                          std::size_t start_j,
                                          std::size_t start_k,
                                          std::size_t end_i, std::size_t end_j,
                                          std::size_t end_k) const override;

 protected:
  std::size_t handleTransformZ(Real z) const override;

  std::size_t handleTransformZWithoutCheck(Real z) const override;
};

}  // namespace xfdtd

#endif  // _XFDTD_CORE_GRID_SPACE_2D_H_
