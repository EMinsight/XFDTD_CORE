#include "grid_space/grid_space_3d.h"

#include <xfdtd/grid_space/grid_space.h>

#include <limits>
#include <memory>
#include <sstream>
#include <utility>
#include <xtensor.hpp>

namespace xfdtd {

GridSpace3D::GridSpace3D(double based_dx, double based_dy, double based_dz,
                         xt::xarray<double> e_node_x,
                         xt::xarray<double> e_node_y,
                         xt::xarray<double> e_node_z)
    : GridSpace{based_dx,
                based_dy,
                based_dz,
                GridSpace::Dimension::THREE,
                std::move(e_node_x),
                std::move(e_node_y),
                std::move(e_node_z)} {}

GridSpace3D::GridSpace3D(
    Type type, GridBox global_box, double based_dx, double based_dy,
    double based_dz, double min_dx, double min_dy, double min_dz,
    xt::xarray<double> e_node_x, xt::xarray<double> e_node_y,
    xt::xarray<double> e_node_z, xt::xarray<double> h_node_x,
    xt::xarray<double> h_node_y, xt::xarray<double> h_node_z,
    xt::xarray<double> e_size_x, xt::xarray<double> e_size_y,
    xt::xarray<double> e_size_z, xt::xarray<double> h_size_x,
    xt::xarray<double> h_size_y, xt::xarray<double> h_size_z)
    : GridSpace{Dimension::THREE,
                type,
                global_box,
                based_dx,
                based_dy,
                based_dz,
                min_dx,
                min_dy,
                min_dz,
                std::move(e_node_x),
                std::move(e_node_y),
                std::move(e_node_z),
                std::move(h_node_x),
                std::move(h_node_y),
                std::move(h_node_z),
                std::move(e_size_x),
                std::move(e_size_y),
                std::move(e_size_z),
                std::move(h_size_x),
                std::move(h_size_y),
                std::move(h_size_z)} {
  _max_x = eNodeX().back();
  _max_y = eNodeY().back();
  _max_z = eNodeZ().back();
  _min_x = eNodeX().front();
  _min_y = eNodeY().front();
  _min_z = eNodeZ().front();
}

std::unique_ptr<GridSpace> GridSpace3D::subGridSpace(
    std::size_t start_i, std::size_t start_j, std::size_t start_k,
    std::size_t end_i, std::size_t end_j, std::size_t end_k) const {
  auto max_i = sizeX();
  auto max_j = sizeY();
  auto max_k = sizeZ();

  if (max_i < start_i || max_j < start_j || max_k < start_k) {
    std::stringstream ss;
    ss << "subGridSpace: start index is out of range: start_i=" << start_i
       << ", start_j=" << start_j << ", start_k=" << start_k;
    throw std::out_of_range(ss.str());
  }

  end_i = std::min(end_i, max_i);
  end_j = std::min(end_j, max_j);
  end_k = std::min(end_k, max_k);

  auto based_dx{basedDx()};
  auto based_dy{basedDy()};
  auto based_dz{basedDz()};
  auto min_dx{minDx()};
  auto min_dy{minDy()};
  auto min_dz{minDz()};
  xt::xarray<double> e_node_x{
      xt::view(eNodeX(), xt::range(start_i, end_i + 1))};
  xt::xarray<double> e_node_y{
      xt::view(eNodeY(), xt::range(start_j, end_j + 1))};
  xt::xarray<double> e_node_z{
      xt::view(eNodeZ(), xt::range(start_k, end_k + 1))};
  xt::xarray<double> h_node_x{xt::view(hNodeX(), xt::range(start_i, end_i))};
  xt::xarray<double> h_node_y{xt::view(hNodeY(), xt::range(start_j, end_j))};
  xt::xarray<double> h_node_z{xt::view(hNodeZ(), xt::range(start_k, end_k))};
  xt::xarray<double> e_size_x{xt::view(eSizeX(), xt::range(start_i, end_i))};
  xt::xarray<double> e_size_y{xt::view(eSizeY(), xt::range(start_j, end_j))};
  xt::xarray<double> e_size_z{xt::view(eSizeZ(), xt::range(start_k, end_k))};
  xt::xarray<double> h_size_x{
      xt::view(hSizeX(), xt::range(start_i, end_i + 1))};
  xt::xarray<double> h_size_y{
      xt::view(hSizeY(), xt::range(start_j, end_j + 1))};
  xt::xarray<double> h_size_z{
      xt::view(hSizeZ(), xt::range(start_k, end_k + 1))};

  auto global_box =
      GridBox{Grid{start_i, start_j, start_k},
              Grid{end_i - start_i, end_j - start_j, end_k - start_k}};

  return std::make_unique<GridSpace3D>(
      type(), global_box, based_dx, based_dy, based_dz, min_dx, min_dy, min_dz,
      std::move(e_node_x), std::move(e_node_y), std::move(e_node_z),
      std::move(h_node_x), std::move(h_node_y), std::move(h_node_z),
      std::move(e_size_x), std::move(e_size_y), std::move(e_size_z),
      std::move(h_size_x), std::move(h_size_y), std::move(h_size_z));
}

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

  _max_x = eNodeX().back();
  _max_y = eNodeY().back();
  _max_z = eNodeZ().back();
  _min_x = eNodeX().front();
  _min_y = eNodeY().front();
  _min_z = eNodeZ().front();
}

}  // namespace xfdtd
