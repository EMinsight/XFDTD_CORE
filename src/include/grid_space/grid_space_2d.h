#ifndef _XFDTD_CORE_GRID_SPACE_2D_H_
#define _XFDTD_CORE_GRID_SPACE_2D_H_

#include <xfdtd/grid_space/grid_space.h>

namespace xfdtd {

class GridSpace2D : public GridSpace {
 public:
  GridSpace2D(GridSpaceRegion region, double based_dx, double based_dy,
              xt::xarray<double> e_node_x, xt::xarray<double> e_node_y);

  GridSpace2D(Type type, GridSpaceRegion region, GridBox global_box,
              double based_dx, double based_dy, double min_dx, double min_dy,
              xt::xarray<double> e_node_x, xt::xarray<double> e_node_y,
              xt::xarray<double> h_node_x, xt::xarray<double> h_node_y,
              xt::xarray<double> e_size_x, xt::xarray<double> e_size_y,
              xt::xarray<double> h_size_x, xt::xarray<double> h_size_y);

  void correctGridSpace() override;

  std::unique_ptr<GridSpace> subGridSpace(std::size_t start_i,
                                          std::size_t start_j,
                                          std::size_t start_k,
                                          std::size_t end_i, std::size_t end_j,
                                          std::size_t end_k) const override;

 protected:
  std::size_t handleTransformZ(double z) const override;

  std::size_t handleTransformZWithoutCheck(double z) const override;
};

}  // namespace xfdtd

#endif  // _XFDTD_CORE_GRID_SPACE_2D_H_
