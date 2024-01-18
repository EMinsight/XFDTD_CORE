#ifndef _XFDTD_LIB_GRID_SPACE_2D_H_
#define _XFDTD_LIB_GRID_SPACE_2D_H_

#include <xfdtd/grid_space/grid_space.h>

namespace xfdtd {

class GridSpace2D : public GridSpace {
 public:
  GridSpace2D(GridSpaceRegion region, double based_dx, double based_dy,
              xt::xarray<double> e_node_x, xt::xarray<double> e_node_y);

  void correctGridSpace() override;

 protected:
  std::size_t handleTransformZ(double z) const override;
};

}  // namespace xfdtd

#endif // _XFDTD_LIB_GRID_SPACE_2D_H_
