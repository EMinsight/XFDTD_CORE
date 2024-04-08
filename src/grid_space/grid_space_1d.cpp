
#include "grid_space/grid_space_1d.h"

#include <xfdtd/common/constant.h>
#include <xfdtd/grid_space/grid_space.h>

namespace xfdtd {

GridSpace1D::GridSpace1D(Real based_dz, Array1D<Real> e_node_z)
    : GridSpace{1,
                1,
                based_dz,
                GridSpace::Dimension::ONE,
                Array1D<Real>{constant::NEG_INF, constant::INF},
                Array1D<Real>{constant::NEG_INF, constant::INF},
                std::move(e_node_z)} {}

void GridSpace1D::correctGridSpace() {
  // correctGridSpaceForOne(basedDz(), eNodeZ(), hNodeZ(), eSizeZ(), hSizeZ());
  // setMinDx(1);
  // setMinDy(1);
  // auto dz_unique{xt::unique(eSizeZ())};
  // auto min_dz{std::numeric_limits<Real>::max()};
  // std::for_each(dz_unique.begin(), dz_unique.end(),
  //               [&min_dz](Real dz) { min_dz = std::min(min_dz, dz); });
  // setMinDz(min_dz);
  throw std::runtime_error("Not implemented");
}

std::size_t GridSpace1D::handleTransformX(Real x) const {
  return (x < std::numeric_limits<Real>::infinity()) ? 0 : 1;
}

std::size_t GridSpace1D::handleTransformY(Real y) const {
  return (y < std::numeric_limits<Real>::infinity()) ? 0 : 1;
}

}  // namespace xfdtd
