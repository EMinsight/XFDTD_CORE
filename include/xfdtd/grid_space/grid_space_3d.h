#ifndef _XFDTD_LIB_GRID_SPACE_3D_H_
#define _XFDTD_LIB_GRID_SPACE_3D_H_

#include <xfdtd/grid_space/grid_space.h>

namespace xfdtd {

class GridSpace3D : public GridSpace {
 public:
  GridSpace3D(GridSpaceRegion region, double based_dx, double based_dy,
              double based_dz, xt::xarray<double> e_node_x,
              xt::xarray<double> e_node_y, xt::xarray<double> e_node_z);

  void correctGridSpace() override;
};

}  // namespace xfdtd

#endif  // _XFDTD_LIB_GRID_SPACE_3D_H_
