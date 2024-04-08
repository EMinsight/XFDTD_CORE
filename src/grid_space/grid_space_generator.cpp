#include <xfdtd/boundary/pml.h>
#include <xfdtd/common/constant.h>
#include <xfdtd/coordinate_system/coordinate_system.h>
#include <xfdtd/grid_space/grid_space_generator.h>
#include <xfdtd/shape/shape.h>

#include <memory>
#include <vector>

#include "grid_space/grid_space_1d.h"
#include "grid_space/grid_space_2d.h"
#include "grid_space/grid_space_3d.h"

namespace xfdtd {

std::unique_ptr<GridSpace> GridSpaceGenerator::generateUniformGridSpace(
    const std::vector<const Shape*>& shapes,
    const std::vector<const Boundary*>& boundaries, Real based_dx,
    Real based_dy, Real based_dz) {
  auto domain = extendDomain(Shape::makeWrappedCube(shapes), boundaries,
                             based_dx, based_dy, based_dz);

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

std::unique_ptr<Cube> GridSpaceGenerator::extendDomain(
    std::unique_ptr<Cube> domain,
    const std::vector<const Boundary*>& boundaries, Real based_dx,
    Real based_dy, Real based_dz) {
  if (boundaries.empty()) {
    return domain;
  }
  auto dimension = decideDimension(domain.get());

  auto min_x = domain->originX();
  auto max_x = domain->endX();
  auto min_y = domain->originY();
  auto max_y = domain->endY();
  auto min_z = domain->originZ();
  auto max_z = domain->endZ();

  for (const auto& b : boundaries) {
    if (auto pml = dynamic_cast<const PML*>(b); pml != nullptr) {
      auto direction = pml->direction();
      auto num = pml->thickness();
      auto main_axis = pml->mainAxis();

      auto space_dim = dimension;
      if (space_dim == GridSpace::Dimension::ONE && main_axis != Axis::XYZ::Z) {
        throw XFDTDGridSpaceException("PML has to be in Z direction");
      }

      if (space_dim == GridSpace::Dimension::TWO && main_axis == Axis::XYZ::Z) {
        throw XFDTDGridSpaceException("PML has to be in X or Y direction");
      }

      if (num < 0) {
        continue;
      }

      switch (direction) {
        case Axis::Direction::XN:
          min_x -= num * based_dx;
          break;
        case Axis::Direction::XP:
          max_x += num * based_dx;
          break;
        case Axis::Direction::YN:
          min_y -= num * based_dy;
          break;
        case Axis::Direction::YP:
          max_y += num * based_dy;
          break;
        case Axis::Direction::ZN:
          min_z -= num * based_dz;
          break;
        case Axis::Direction::ZP:
          max_z += num * based_dz;
          break;
        default:
          throw XFDTDGridSpaceException{"Invalid direction"};
      }
    }
  }

  return std::make_unique<Cube>(
      Vector{min_x, min_y, min_z},
      Vector{max_x - min_x, max_y - min_y, max_z - min_z});
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
    const Cube* domain, Real dz) {
  std::size_t nx{1};
  std::size_t ny{1};
  std::size_t nz{1};
  auto min_z{std::numeric_limits<Real>::max()};
  auto max_z{std::numeric_limits<Real>::min()};

  auto cube{domain->wrappedCube()};
  min_z = std::min(min_z, cube->originZ());
  max_z = std::max(max_z, cube->endZ());

  nz = std::round<size_t>((max_z - min_z) / dz);
  if (nz <= 1) {
    throw XFDTDGridSpaceException{"GridSpace is too small"};
  }

  max_z = min_z + dz * nz;

  auto region{Cube{{constant::NEG_INF, constant::NEG_INF, min_z},
                   {constant::INF, constant::INF, nz * dz}}};
  auto e_node_z{xt::linspace<Real>(region.originZ(), region.endZ(), nz + 1)};

  return std::make_unique<GridSpace1D>(dz, std::move(e_node_z));
}

std::unique_ptr<GridSpace> GridSpaceGenerator::generateGridSpace2D(
    const Cube* domain, Real dx, Real dy) {
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

  max_x = min_x + dx * nx;
  max_y = min_y + dy * ny;

  auto region{Cube{{min_x, min_y, constant::NEG_INF},
                   {nx * dx, ny * dy, constant::INF}}};

  auto e_node_x{xt::linspace<Real>(region.originX(), region.endX(), nx + 1)};
  auto e_node_y{xt::linspace<Real>(region.originY(), region.endY(), ny + 1)};

  return std::make_unique<GridSpace2D>(dx, dy, std::move(e_node_x),
                                       std::move(e_node_y));
}

std::unique_ptr<GridSpace> GridSpaceGenerator::generateGridSpace3D(
    const Cube* domain, Real dx, Real dy, Real dz) {
  std::size_t nx{1};
  std::size_t ny{1};
  std::size_t nz{1};

  auto min_x = domain->originX();
  auto max_x = domain->endX();
  auto min_y = domain->originY();
  auto max_y = domain->endY();
  auto min_z = domain->originZ();
  auto max_z = domain->endZ();

  nx = std::round<size_t>((max_x - min_x) / dx);
  ny = std::round<size_t>((max_y - min_y) / dy);
  nz = std::round<size_t>((max_z - min_z) / dz);
  if (nx <= 1 || ny <= 1 || nz <= 1) {
    throw XFDTDGridSpaceException{"GridSpace is too small"};
  }

  max_x = min_x + dx * nx;
  max_y = min_y + dy * ny;
  max_z = min_z + dz * nz;

  auto region{Cube{{min_x, min_y, min_z},
                   {max_x - min_x, max_y - min_y, max_z - min_z}}};

  auto e_node_x{xt::linspace<Real>(region.originX(), region.endX(), nx + 1)};
  auto e_node_y{xt::linspace<Real>(region.originY(), region.endY(), ny + 1)};
  auto e_node_z{xt::linspace<Real>(region.originZ(), region.endZ(), nz + 1)};

  return std::make_unique<GridSpace3D>(dx, dy, dz, std::move(e_node_x),
                                       std::move(e_node_y),
                                       std::move(e_node_z));
}

}  // namespace xfdtd