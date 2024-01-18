#include <xfdtd/grid_space/grid_space_3d.h>

#include <limits>
#include <utility>
#include <xtensor.hpp>

namespace xfdtd {

GridSpace3D::GridSpace3D(GridSpaceRegion region, double based_dx,
                         double based_dy, double based_dz,
                         xt::xarray<double> e_node_x,
                         xt::xarray<double> e_node_y,
                         xt::xarray<double> e_node_z)
    : GridSpace{std::move(region),
                based_dx,
                based_dy,
                based_dz,
                GridSpace::Dimension::THREE,
                std::move(e_node_x),
                std::move(e_node_y),
                std::move(e_node_z)} {}

void GridSpace3D::correctGridSpace() {
  correctGridSpaceForOne(basedDx(), eNodeX(), hNodeX(), eSizeX(), hSizeX());
  correctGridSpaceForOne(basedDy(), eNodeY(), hNodeY(), eSizeY(), hSizeY());
  correctGridSpaceForOne(basedDz(), eNodeZ(), hNodeZ(), eSizeZ(), hSizeZ());
  auto dx_unique{xt::unique(eSizeX())};
  auto min_dx{std::numeric_limits<double>::max()};
  std::for_each(dx_unique.begin(), dx_unique.end(),
                [&min_dx](double dx) { min_dx = std::min(min_dx, dx); });
  setMinDx(min_dx);
  auto dy_unique{xt::unique(eSizeY())};
  auto min_dy{std::numeric_limits<double>::max()};
  std::for_each(dy_unique.begin(), dy_unique.end(),
                [&min_dy](double dy) { min_dy = std::min(min_dy, dy); });
  setMinDy(min_dy);
  auto dz_unique{xt::unique(eSizeZ())};
  auto min_dz{std::numeric_limits<double>::max()};
  std::for_each(dz_unique.begin(), dz_unique.end(),
                [&min_dz](double dz) { min_dz = std::min(min_dz, dz); });
  setMinDz(min_dz);

  generateGrid(hNodeX().size(), hNodeY().size(), hNodeZ().size());
  _nx = hNodeX().size();
  _ny = hNodeY().size();
  _nz = hNodeZ().size();
  _max_x = eNodeX().back();
  _max_y = eNodeY().back();
  _max_z = eNodeZ().back();
  _min_x = eNodeX().front();
  _min_y = eNodeY().front();
  _min_z = eNodeZ().front();
}

}  // namespace xfdtd
