#include "grid_space/grid_space_2d.h"

#include <limits>
#include <memory>
#include <utility>
#include <xtensor.hpp>

#include "xfdtd/grid_space/grid_space.h"
#include "xfdtd/util/constant.h"

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

GridSpace2D::GridSpace2D(
    Type type, GridSpaceRegion region, GridBox global_box, double based_dx,
    double based_dy, double min_dx, double min_dy, xt::xarray<double> e_node_x,
    xt::xarray<double> e_node_y, xt::xarray<double> h_node_x,
    xt::xarray<double> h_node_y, xt::xarray<double> e_size_x,
    xt::xarray<double> e_size_y, xt::xarray<double> h_size_x,
    xt::xarray<double> h_size_y)
    : GridSpace(
          GridSpace::Dimension::TWO, type, std::move(region), global_box,
          based_dx, based_dy, 1, min_dx, min_dy, 1, std::move(e_node_x),
          std::move(e_node_y),
          xt::xarray<double>({constant::NEG_INF, constant::INF}),
          std::move(h_node_x), std::move(h_node_y), xt::xarray<double>({0}),
          std::move(e_size_x), std::move(e_size_y), xt::xarray<double>({1}),
          std::move(h_size_x), std::move(h_size_y), xt::xarray<double>({1})) {
  _max_x = eNodeX().back();
  _max_y = eNodeY().back();
  _max_z = constant::INF;
  _min_x = eNodeX().front();
  _min_y = eNodeY().front();
  _min_z = constant::NEG_INF;
}

std::unique_ptr<GridSpace> GridSpace2D::subGridSpace(
    std::size_t start_i, std::size_t start_j, std::size_t start_k,
    std::size_t end_i, std::size_t end_j, std::size_t end_k) const {
  auto max_i = sizeX();
  auto max_j = sizeY();
  auto max_k = sizeZ();

  if (max_i < start_i || max_j < start_j || max_k < start_k) {
    throw XFDTDGridSpaceException("subGridSpace: start index is out of range");
  }

  end_i = std::min(end_i, max_i);
  end_j = std::min(end_j, max_j);
  end_k = std::min(end_k, max_k);

  auto based_dx{basedDx()};
  auto based_dy{basedDy()};
  auto min_dx{minDx()};
  auto min_dy{minDy()};
  xt::xarray<double> e_node_x{
      xt::view(eNodeX(), xt::range(start_i, end_i + 1))};
  xt::xarray<double> e_node_y{
      xt::view(eNodeY(), xt::range(start_j, end_j + 1))};

  xt::xarray<double> h_node_x{xt::view(hNodeX(), xt::range(start_i, end_i))};
  xt::xarray<double> h_node_y{xt::view(hNodeY(), xt::range(start_j, end_j))};
  xt::xarray<double> e_size_x{xt::view(eSizeX(), xt::range(start_i, end_i))};
  xt::xarray<double> e_size_y{xt::view(eSizeY(), xt::range(start_j, end_j))};
  xt::xarray<double> h_size_x{
      xt::view(hSizeX(), xt::range(start_i, end_i + 1))};
  xt::xarray<double> h_size_y{
      xt::view(hSizeY(), xt::range(start_j, end_j + 1))};

  auto region = GridSpaceRegion{
      Vector{e_node_x.front(), e_node_y.front(), constant::NEG_INF},
      Vector{e_node_x.back(), e_node_y.back(), constant::INF}};

  auto global_box = GridBox{Grid{start_i, start_j, 0},
                            Grid{end_i - start_i, end_j - start_j, 1}};

  return std::make_unique<GridSpace2D>(
      type(), std::move(region), global_box, based_dx, based_dy, min_dx, min_dy,
      std::move(e_node_x), std::move(e_node_y), std::move(h_node_x),
      std::move(h_node_y), std::move(e_size_x), std::move(e_size_y),
      std::move(h_size_x), std::move(h_size_y));
}

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
  setMinDz(1.0);
  eSizeZ() = xt::zeros<double>({1});
  eSizeZ()(0) = 1;
  hSizeZ() = xt::zeros<double>({1});
  hSizeZ()(0) = 1;
  hNodeZ() = xt::zeros<double>({1});

  _max_x = eNodeX().back();
  _max_y = eNodeY().back();
  _max_z = constant::INF;
  _min_x = eNodeX().front();
  _min_y = eNodeY().front();
  _min_z = constant::NEG_INF;
}

std::size_t GridSpace2D::handleTransformZWithoutCheck(double z) const {
  return (z < std::numeric_limits<double>::infinity()) ? 0 : 1;
}

std::size_t GridSpace2D::handleTransformZ(double z) const {
  return (z < std::numeric_limits<double>::infinity()) ? 0 : 1;
}

}  // namespace xfdtd
