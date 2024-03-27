#ifndef _XFDTD_CORE_GRID_SPACE_3D_H_
#define _XFDTD_CORE_GRID_SPACE_3D_H_

#include <xfdtd/grid_space/grid_space.h>

namespace xfdtd {

class GridSpace3D : public GridSpace {
 public:
  GridSpace3D(double based_dx, double based_dy, double based_dz,
              xt::xarray<double> e_node_x, xt::xarray<double> e_node_y,
              xt::xarray<double> e_node_z);

  GridSpace3D(Type type, GridBox global_box, double based_dx, double based_dy,
              double based_dz, double min_dx, double min_dy, double min_dz,
              xt::xarray<double> e_node_x, xt::xarray<double> e_node_y,
              xt::xarray<double> e_node_z, xt::xarray<double> h_node_x,
              xt::xarray<double> h_node_y, xt::xarray<double> h_node_z,
              xt::xarray<double> e_size_x, xt::xarray<double> e_size_y,
              xt::xarray<double> e_size_z, xt::xarray<double> h_size_x,
              xt::xarray<double> h_size_y, xt::xarray<double> h_size_z);

  std::unique_ptr<GridSpace> subGridSpace(std::size_t start_i,
                                          std::size_t start_j,
                                          std::size_t start_k,
                                          std::size_t end_i, std::size_t end_j,
                                          std::size_t end_k) const override;

  void correctGridSpace() override;
};

}  // namespace xfdtd

#endif  // _XFDTD_CORE_GRID_SPACE_3D_H_
