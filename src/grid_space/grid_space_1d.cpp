#include "grid_space/grid_space_1d.h"

#include <xfdtd/common/constant.h>
#include <xfdtd/common/type_define.h>
#include <xfdtd/grid_space/grid_space.h>

#include <sstream>
#include <xtensor.hpp>

namespace xfdtd {

GridSpace1D::GridSpace1D(Real based_dz, Array1D<Real> e_node_z)
    : GridSpace{1,
                1,
                based_dz,
                GridSpace::Dimension::ONE,
                Array1D<Real>{constant::NEG_INF, constant::INF},
                Array1D<Real>{constant::NEG_INF, constant::INF},
                std::move(e_node_z)} {}

GridSpace1D::GridSpace1D(Type type, GridBox global_box, Real based_dz,
                         Real min_dz, Array1D<Real> e_node_z,
                         Array1D<Real> h_node_z, Array1D<Real> e_size_z,
                         Array1D<Real> h_size_z)
    : GridSpace{Dimension::ONE,
                type,
                global_box,
                1,
                1,
                based_dz,
                1,
                1,
                min_dz,
                Array1D<Real>{constant::NEG_INF, constant::INF},
                Array1D<Real>{constant::NEG_INF, constant::INF},
                std::move(e_node_z),
                Array1D<Real>{0},
                Array1D<Real>{0},
                std::move(h_node_z),
                Array1D<Real>{1},
                Array1D<Real>{1},
                std::move(e_size_z),
                Array1D<Real>{1},
                Array1D<Real>{1},
                std::move(h_size_z)} {
  _max_x = constant::INF;
  _max_y = constant::INF;
  _max_z = eNodeZ().back();
  _min_x = constant::NEG_INF;
  _min_y = constant::NEG_INF;
  _min_z = eNodeZ().front();
}

void GridSpace1D::correctGridSpace() {
  correctGridSpaceForOne(basedDz(), eNodeZ(), hNodeZ(), eSizeZ(), hSizeZ());
  setMinDx(1);
  setMinDy(1);
  auto dz_unique{xt::unique(eSizeZ())};
  auto min_dz{std::numeric_limits<Real>::max()};
  std::for_each(dz_unique.begin(), dz_unique.end(),
                [&min_dz](Real dz) { min_dz = std::min(min_dz, dz); });
  setMinDz(min_dz);

  eSizeX() = xt::zeros<Real>({1});
  eSizeX()(0) = 1;
  hSizeX() = xt::zeros<Real>({1});
  hSizeX()(0) = 1;
  hNodeX() = xt::zeros<Real>({1});

  eSizeY() = xt::zeros<Real>({1});
  eSizeY()(0) = 1;
  hSizeY() = xt::zeros<Real>({1});
  hSizeY()(0) = 1;
  hNodeY() = xt::zeros<Real>({1});

  _max_x = constant::INF;
  _max_y = constant::INF;
  _max_z = eNodeZ().back();
  _min_x = constant::NEG_INF;
  _min_y = constant::NEG_INF;
  _min_z = eNodeZ().front();
}

auto GridSpace1D::subGridSpace(Index start_i, Index start_j, Index start_k,
                               Index end_i, Index end_j, Index end_k) const
    -> std::unique_ptr<GridSpace> {
  auto max_i = sizeX();
  auto max_j = sizeY();
  auto max_k = sizeZ();

  if (max_i < start_i || max_j < start_j || max_k < start_k) {
    std::stringstream ss;
    ss << "Invalid sub grid space: start_i=" << start_i
       << ", start_j=" << start_j << ", start_k=" << start_k
       << ", max_i=" << max_i << ", max_j=" << max_j << ", max_k=" << max_k;
    throw XFDTDGridSpaceException(ss.str());
  }

  end_i = std::min(end_i, max_i);
  end_j = std::min(end_j, max_j);
  end_k = std::min(end_k, max_k);

  auto based_dz{basedDz()};
  auto min_dz{minDz()};
  Array1D<Real> e_node_z{xt::view(eNodeZ(), xt::range(start_k, end_k + 1))};
  Array1D<Real> h_node_z{xt::view(hNodeZ(), xt::range(start_k, end_k))};
  Array1D<Real> e_size_z{xt::view(eSizeZ(), xt::range(start_k, end_k))};
  Array1D<Real> h_size_z{xt::view(hSizeZ(), xt::range(start_k, end_k + 1))};

  auto global_box =
      GridBox{Grid{start_i, start_j, start_k},
              Grid{end_i - start_i, end_j - start_j, end_k - start_k}};
  return std::make_unique<GridSpace1D>(
      type(), global_box, based_dz, min_dz, std::move(e_node_z),
      std::move(h_node_z), std::move(e_size_z), std::move(h_size_z));
};

static inline auto compareWithINF(Real x) -> bool { return x < constant::INF; }

auto GridSpace1D::handleTransformXWithoutCheck(Real x) const -> std::size_t {
  return compareWithINF(x) ? 0 : 1;
}

auto GridSpace1D::handleTransformYWithoutCheck(Real y) const -> std::size_t {
  return compareWithINF(y) ? 0 : 1;
}

std::size_t GridSpace1D::handleTransformX(Real x) const {
  return compareWithINF(x) ? 0 : 1;
}

std::size_t GridSpace1D::handleTransformY(Real y) const {
  return compareWithINF(y) ? 0 : 1;
}

}  // namespace xfdtd
