#ifndef _XFDTD_CORE_GRID_SPACE_1D_H_
#define _XFDTD_CORE_GRID_SPACE_1D_H_

#include <xfdtd/grid_space/grid_space.h>

#include <stdexcept>

namespace xfdtd {

class GridSpace1D : public GridSpace {
 public:
  GridSpace1D(GridSpaceRegion region, double based_dz,
              xt::xarray<double> e_node_z);

  void correctGridSpace() override;

  std::unique_ptr<GridSpace> subGridSpace(std::size_t start_i,
                                          std::size_t start_j,
                                          std::size_t start_k,
                                          std::size_t end_i, std::size_t end_j,
                                          std::size_t end_k) const override {
    throw std::runtime_error("subGridSpace is not implemented");
  }

 protected:
  std::size_t handleTransformXWithoutCheck(double x) const override {
    throw std::runtime_error("handleTransformXWithoutCheck is not implemented");
  }

  std::size_t handleTransformYWithoutCheck(double y) const override {
    throw std::runtime_error("handleTransformYWithoutCheck is not implemented");
  }

  std::size_t handleTransformX(double x) const override;

  std::size_t handleTransformY(double y) const override;
};

}  // namespace xfdtd

#endif  // _XFDTD_CORE_GRID_SPACE_1D_H_
