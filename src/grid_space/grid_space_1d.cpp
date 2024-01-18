
#include <xfdtd/grid_space/grid_space.h>
#include <xfdtd/grid_space/grid_space_1d.h>

#include <limits>
#include <xtensor.hpp>

namespace xfdtd {

GridSpace1D::GridSpace1D(GridSpaceRegion region, double based_dz,
                         xt::xarray<double> e_node_z)
    : GridSpace{std::move(region),
                1,
                1,
                based_dz,
                GridSpace::Dimension::ONE,
                xt::xarray<double>{-1 * std::numeric_limits<double>::infinity(),
                                   std::numeric_limits<double>::infinity()},
                xt::xarray<double>{-1 * std::numeric_limits<double>::infinity(),
                                   std::numeric_limits<double>::infinity()},
                std::move(e_node_z)} {}

void GridSpace1D::correctGridSpace() {
  correctGridSpaceForOne(basedDz(), eNodeZ(), hNodeZ(), eSizeZ(), hSizeZ());
  setMinDx(1);
  setMinDy(1);
  auto dz_unique{xt::unique(eSizeZ())};
  auto min_dz{std::numeric_limits<double>::max()};
  std::for_each(dz_unique.begin(), dz_unique.end(),
                [&min_dz](double dz) { min_dz = std::min(min_dz, dz); });
  setMinDz(min_dz);
}

std::size_t GridSpace1D::handleTransformX(double x) const {
  return (x < std::numeric_limits<double>::infinity()) ? 0 : 1;
}

std::size_t GridSpace1D::handleTransformY(double y) const {
  return (y < std::numeric_limits<double>::infinity()) ? 0 : 1;
}

}  // namespace xfdtd
