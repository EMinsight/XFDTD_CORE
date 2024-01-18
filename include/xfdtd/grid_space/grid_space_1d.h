#ifndef _XFDTD_LIB_GRID_SPACE_1D_H_
#define _XFDTD_LIB_GRID_SPACE_1D_H_

#include <xfdtd/grid_space/grid_space.h>

namespace xfdtd {

class GridSpace1D : public GridSpace {
 public:
  GridSpace1D(GridSpaceRegion region, double based_dz,
              xt::xarray<double> e_node_z);

  void correctGridSpace() override;

 protected:
  std::size_t handleTransformX(double x) const override;

  std::size_t handleTransformY(double y) const override;
};

}  // namespace xfdtd

#endif // _XFDTD_LIB_GRID_SPACE_1D_H_
