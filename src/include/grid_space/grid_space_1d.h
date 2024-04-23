#ifndef _XFDTD_CORE_GRID_SPACE_1D_H_
#define _XFDTD_CORE_GRID_SPACE_1D_H_

#include <xfdtd/grid_space/grid_space.h>

namespace xfdtd {

class GridSpace1D : public GridSpace {
 public:
  GridSpace1D(Real based_dz, Array1D<Real> e_node_z);

  GridSpace1D(Type type, GridBox global_box, Real based_dz, Real min_dz,
              Array1D<Real> e_node_z, Array1D<Real> h_node_z,
              Array1D<Real> e_size_z, Array1D<Real> h_size_z);

  void correctGridSpace() override;

  std::unique_ptr<GridSpace> subGridSpace(std::size_t start_i,
                                          std::size_t start_j,
                                          std::size_t start_k,
                                          std::size_t end_i, std::size_t end_j,
                                          std::size_t end_k) const override;

  ~GridSpace1D() override = default;

 protected:
  std::size_t handleTransformXWithoutCheck(Real x) const override;

  std::size_t handleTransformYWithoutCheck(Real y) const override;

  std::size_t handleTransformX(Real x) const override;

  std::size_t handleTransformY(Real y) const override;
};

}  // namespace xfdtd

#endif  // _XFDTD_CORE_GRID_SPACE_1D_H_
