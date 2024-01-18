#include <xfdtd/grid_space/grid_space_2d.h>

#include <limits>
#include <utility>
#include <xtensor.hpp>

namespace xfdtd {

GridSpace2D::GridSpace2D(GridSpaceRegion region, double based_dx,
                         double based_dy, xt::xarray<double> e_node_x,
                         xt::xarray<double> e_node_y)
    : GridSpace{std::move(region),
                based_dx,
                based_dx,
                1.0,
                GridSpace::Dimension::TWO,
                std::move(e_node_x),
                std::move(e_node_y),
                {-1 * std::numeric_limits<double>::infinity(),
                 std::numeric_limits<double>::infinity()}} {}

void GridSpace2D::correctGridSpace() {
  correctGridSpaceForOne(basedDx(), eNodeX(), hNodeX(), eSizeX(), hSizeX());
  correctGridSpaceForOne(basedDy(), eNodeY(), hNodeY(), eSizeY(), hSizeY());
  setMinDz(1);
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
  eSizeZ() = xt::zeros<double>({1});
  eSizeZ()(0) = 1;
  hSizeZ() = xt::zeros<double>({1});
  hSizeZ()(0) = 1;
  hNodeZ() = xt::zeros<double>({1});
  generateGrid(hNodeX().size(), hNodeY().size(), hNodeZ().size());
}

std::size_t GridSpace2D::handleTransformZ(double z) const {
  return (z < std::numeric_limits<double>::infinity()) ? 0 : 1;
}

}  // namespace xfdtd
