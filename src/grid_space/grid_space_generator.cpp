
#include <xfdtd/boundary/pml.h>
#include <xfdtd/grid_space/grid_space_generator.h>
#include <xfdtd/shape/shape.h>

#include <cmath>
#include <limits>
#include <memory>
#include <vector>
#include <xtensor.hpp>

#include "grid_space/grid_space_1d.h"
#include "grid_space/grid_space_2d.h"
#include "grid_space/grid_space_3d.h"

namespace xfdtd {

std::unique_ptr<GridSpace> GridSpaceGenerator::generateUniformGridSpace(
    const std::vector<const Shape*>& shapes, double based_dx, double based_dy,
    double based_dz) {
  auto domain = Shape::makeWrappedCube(shapes);

  auto dimension{decideDimension(domain.get())};

  switch (dimension) {
    case GridSpace::Dimension::ONE:
      return generateGridSpace1D(domain.get(), based_dz);
    case GridSpace::Dimension::TWO:
      return generateGridSpace2D(domain.get(), based_dx, based_dy);
    case GridSpace::Dimension::THREE:
      return generateGridSpace3D(domain.get(), based_dx, based_dy, based_dz);
    default:
      throw XFDTDGridSpaceException{"GridSpace dimension is undefined"};
  }
}

GridSpace::Dimension GridSpaceGenerator::decideDimension(const Shape* shape) {
  auto dimension{GridSpace::Dimension::UNDEFINED};

  auto cube{shape->wrappedCube()};
  if (std::isinf(cube->sizeX()) && std::isinf(cube->sizeY()) &&
      !std::isinf(cube->sizeZ())) {
    dimension = GridSpace::Dimension::ONE;
    return dimension;
  }

  if (!std::isinf(cube->sizeX()) && !std::isinf(cube->sizeY()) &&
      std::isinf(cube->sizeZ())) {
    dimension = GridSpace::Dimension::TWO;
    return dimension;
  }

  if (!std::isinf(cube->sizeX()) && !std::isinf(cube->sizeY()) &&
      !std::isinf(cube->sizeZ())) {
    dimension = GridSpace::Dimension::THREE;
    return dimension;
  }

  return dimension;
}

std::unique_ptr<GridSpace> GridSpaceGenerator::generateGridSpace1D(
    const Cube* domain, double dz) {
  std::size_t nx{1};
  std::size_t ny{1};
  std::size_t nz{1};
  auto min_z{std::numeric_limits<double>::max()};
  auto max_z{std::numeric_limits<double>::min()};

  auto cube{domain->wrappedCube()};
  min_z = std::min(min_z, cube->originZ());
  max_z = std::max(max_z, cube->endZ());

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
    const Cube* domain, double dx, double dy) {
  std::size_t nx{1};
  std::size_t ny{1};
  std::size_t nz{1};

  auto min_x = domain->originX();
  auto max_x = domain->endX();
  auto min_y = domain->originY();
  auto max_y = domain->endY();

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
    const Cube* domain, double dx, double dy, double dz) {
  std::size_t nx{1};
  std::size_t ny{1};
  std::size_t nz{1};

  auto min_x{std::numeric_limits<double>::max()};
  auto max_x{std::numeric_limits<double>::min()};
  auto min_y{std::numeric_limits<double>::max()};
  auto max_y{std::numeric_limits<double>::min()};
  auto min_z{std::numeric_limits<double>::max()};
  auto max_z{std::numeric_limits<double>::min()};

  min_x = domain->originX();
  max_x = domain->endX();
  min_y = domain->originY();
  max_y = domain->endY();
  min_z = domain->originZ();
  max_z = domain->endZ();

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