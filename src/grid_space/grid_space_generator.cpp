
#include "grid_space/grid_space_generator.h"

#include <xfdtd/grid_space/grid_space.h>
#include <xfdtd/grid_space/grid_space_1d.h>
#include <xfdtd/grid_space/grid_space_2d.h>
#include <xfdtd/grid_space/grid_space_3d.h>

#include <cmath>
#include <limits>
#include <memory>
#include <xtensor.hpp>

#include "xfdtd/boundary/pml.h"

namespace xfdtd {

std::unique_ptr<GridSpace> GridSpaceGenerator::generate(
    const std::vector<std::shared_ptr<Shape>>& shapes, double based_dx,
    double based_dy, double based_dz) {
  auto space_dimension{decideDimension(shapes)};

  if (space_dimension == GridSpace::Dimension::UNDEFINED) {
    throw XFDTDGridSpaceException{"GridSpace dimension is undefined"};
  }

  std::unique_ptr<GridSpace> g{nullptr};

  switch (space_dimension) {
    case GridSpace::Dimension::ONE:
      g = generateGridSpace1D(shapes, based_dz);
      break;
    case GridSpace::Dimension::TWO:
      g = generateGridSpace2D(shapes, based_dx, based_dy);
      break;
    case GridSpace::Dimension::THREE:
      g = generateGridSpace3D(shapes, based_dx, based_dy, based_dz);
      break;
    default:
      throw XFDTDGridSpaceException{"GridSpace dimension is undefined"};
  }

  if (g == nullptr) {
    throw XFDTDGridSpaceException{"GridSpace is not generated"};
  }

  g->correctGridSpace();

  return g;
}

std::unique_ptr<GridSpace> GridSpaceGenerator::generate(
    const std::vector<std::shared_ptr<Shape>>& shapes,
    const std::vector<std::shared_ptr<Boundary>>& boundaries, double based_dx,
    double based_dy, double based_dz) {
  auto space_dimension{decideDimension(shapes)};

  if (space_dimension == GridSpace::Dimension::UNDEFINED) {
    throw XFDTDGridSpaceException{"GridSpace dimension is undefined"};
  }

  std::unique_ptr<GridSpace> g{nullptr};

  switch (space_dimension) {
    case GridSpace::Dimension::ONE:
      g = generateGridSpace1D(shapes, based_dz);
      break;
    case GridSpace::Dimension::TWO:
      g = generateGridSpace2D(shapes, based_dx, based_dy);
      break;
    case GridSpace::Dimension::THREE:
      g = generateGridSpace3D(shapes, based_dx, based_dy, based_dz);
      break;
    default:
      throw XFDTDGridSpaceException{"GridSpace dimension is undefined"};
  }
  g->correctGridSpace(); // TODO(franzero): 

  if (g == nullptr) {
    throw XFDTDGridSpaceException{"GridSpace is not generated"};
  }

  for (const auto& b : boundaries) {
    auto pml{std::dynamic_pointer_cast<PML>(b)};
    if (pml == nullptr) {
      throw XFDTDGridSpaceException{"Only PML is supported for now"};
    }

    auto thickness{pml->thickness()};
    if (thickness <= 0) {
      continue;
    }

    auto direction{pml->direction()};
    auto main_axis{pml->mainAxis()};
    switch (main_axis) {
      case Axis::XYZ::X:
        g->extendGridSpace(direction, thickness, based_dx);
        continue;
      case Axis::XYZ::Y:
        g->extendGridSpace(direction, thickness, based_dy);
        continue;
      case Axis::XYZ::Z:
        g->extendGridSpace(direction, thickness, based_dz);
        continue;
      default:
        throw XFDTDGridSpaceException{"Invalid main axis"};
    }
  }

  g->correctGridSpace();

  return g;
}

GridSpace::Dimension GridSpaceGenerator::decideDimension(
    const std::vector<std::shared_ptr<Shape>>& shapes) {
  auto dimension{GridSpace::Dimension::UNDEFINED};
  for (const auto& shape : shapes) {
    auto cube{shape->wrappedCube()};
    if (std::isinf(cube.sizeX()) && std::isinf(cube.sizeY()) &&
        !std::isinf(cube.sizeZ())) {
      if (dimension == GridSpace::Dimension::TWO ||
          dimension == GridSpace::Dimension::THREE) {
        throw XFDTDGridSpaceException{
            "Dimension One is wrong for having TWO or THREE dimension shapes"};
      }
      dimension = GridSpace::Dimension::ONE;
    }

    if (!std::isinf(cube.sizeX()) && !std::isinf(cube.sizeY()) &&
        std::isinf(cube.sizeZ())) {
      if (dimension == GridSpace::Dimension::THREE) {
        throw XFDTDGridSpaceException{
            "Dimension Two is wrong for having THREE dimension shapes"};
      }
      dimension = GridSpace::Dimension::TWO;
    }

    if (!std::isinf(cube.sizeX()) && !std::isinf(cube.sizeY()) &&
        !std::isinf(cube.sizeZ())) {
      dimension = GridSpace::Dimension::THREE;
    }
  }
  return dimension;
}

std::unique_ptr<GridSpace> GridSpaceGenerator::generateGridSpace1D(
    const std::vector<std::shared_ptr<Shape>>& shapes, double dz) {
  std::size_t nx{1};
  std::size_t ny{1};
  std::size_t nz{1};
  auto min_z{std::numeric_limits<double>::max()};
  auto max_z{std::numeric_limits<double>::min()};
  for (const auto& shape : shapes) {
    auto cube{shape->wrappedCube()};
    min_z = std::min(min_z, cube.originZ());
    max_z = std::max(max_z, cube.endZ());
  }
  nz = std::round<size_t>((max_z - min_z) / dz);
  if (nz <= 1) {
    throw XFDTDGridSpaceException{"GridSpace is too small"};
  }

  auto region{GridSpace::GridSpaceRegion{
      {-1 * std::numeric_limits<double>::infinity(),
       -1 * std::numeric_limits<double>::infinity(), min_z},
      {std::numeric_limits<double>::infinity(),
       std::numeric_limits<double>::infinity(), nz * dz}}};
  auto e_node_z{xt::linspace<double>(region.originZ(), region.endZ(), nz + 1)};

  return std::make_unique<GridSpace1D>(std::move(region), dz,
                                       std::move(e_node_z));
}

std::unique_ptr<GridSpace> GridSpaceGenerator::generateGridSpace2D(
    const std::vector<std::shared_ptr<Shape>>& shapes, double dx, double dy) {
  std::size_t nx{1};
  std::size_t ny{1};
  std::size_t nz{1};

  auto min_x{std::numeric_limits<double>::max()};
  auto max_x{std::numeric_limits<double>::min()};
  auto min_y{std::numeric_limits<double>::max()};
  auto max_y{std::numeric_limits<double>::min()};

  for (const auto& shape : shapes) {
    auto cube{shape->wrappedCube()};
    min_x = std::min(min_x, cube.originX());
    max_x = std::max(max_x, cube.endX());
    min_y = std::min(min_y, cube.originY());
    max_y = std::max(max_y, cube.endY());
  }

  nx = std::round<size_t>((max_x - min_x) / dx);
  ny = std::round<size_t>((max_y - min_y) / dy);
  if (nx <= 1 || ny <= 1) {
    throw XFDTDGridSpaceException{"GridSpace is too small"};
  }

  auto region{GridSpace::GridSpaceRegion{
      {min_x, min_y, -1 * std::numeric_limits<double>::infinity()},
      {nx * dx, ny * dy, std::numeric_limits<double>::infinity()}}};

  auto e_node_x{xt::linspace<double>(region.originX(), region.endX(), nx + 1)};
  auto e_node_y{xt::linspace<double>(region.originY(), region.endY(), ny + 1)};

  return std::make_unique<GridSpace2D>(
      std::move(region), dx, dy, std::move(e_node_x), std::move(e_node_y));
}

std::unique_ptr<GridSpace> GridSpaceGenerator::generateGridSpace3D(
    const std::vector<std::shared_ptr<Shape>>& shapes, double dx, double dy,
    double dz) {
  std::size_t nx{1};
  std::size_t ny{1};
  std::size_t nz{1};

  auto min_x{std::numeric_limits<double>::max()};
  auto max_x{std::numeric_limits<double>::min()};
  auto min_y{std::numeric_limits<double>::max()};
  auto max_y{std::numeric_limits<double>::min()};
  auto min_z{std::numeric_limits<double>::max()};
  auto max_z{std::numeric_limits<double>::min()};

  for (const auto& shape : shapes) {
    auto cube{shape->wrappedCube()};
    min_x = std::min(min_x, cube.originX());
    max_x = std::max(max_x, cube.endX());
    min_y = std::min(min_y, cube.originY());
    max_y = std::max(max_y, cube.endY());
    min_z = std::min(min_z, cube.originZ());
    max_z = std::max(max_z, cube.endZ());
  }

  nx = std::round<size_t>((max_x - min_x) / dx);
  ny = std::round<size_t>((max_y - min_y) / dy);
  nz = std::round<size_t>((max_z - min_z) / dz);
  if (nx <= 1 || ny <= 1 || nz <= 1) {
    throw XFDTDGridSpaceException{"GridSpace is too small"};
  }

  auto region{GridSpace::GridSpaceRegion{{min_x, min_y, min_z},
                                         {nx * dx, ny * dy, nz * dz}}};

  auto e_node_x{xt::linspace<double>(region.originX(), region.endX(), nx + 1)};
  auto e_node_y{xt::linspace<double>(region.originY(), region.endY(), ny + 1)};
  auto e_node_z{xt::linspace<double>(region.originZ(), region.endZ(), nz + 1)};

  return std::make_unique<GridSpace3D>(std::move(region), dx, dy, dz,
                                       std::move(e_node_x), std::move(e_node_y),
                                       std::move(e_node_z));
}

}  // namespace xfdtd