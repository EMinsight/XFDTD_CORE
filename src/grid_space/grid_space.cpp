#include <xfdtd/grid_space/grid_space.h>

#include <memory>
#include <utility>
#include <xtensor.hpp>

#include "util/float_compare.h"
#include "xfdtd/coordinate_system/coordinate_system.h"
#include "xfdtd/shape/shape.h"

namespace xfdtd {

Grid::Grid(std::size_t i, std::size_t j, std::size_t k) : _i{i}, _j{j}, _k{k} {}

Grid Grid::operator+(const Grid& grid) const {
  return {_i + grid._i, _j + grid._j, _k + grid._k};
}

Grid Grid::operator-(const Grid& grid) const {
  return {_i - grid._i, _j - grid._j, _k - grid._k};
}

std::size_t Grid::i() const { return _i; }

std::size_t Grid::j() const { return _j; }

std::size_t Grid::k() const { return _k; }

GridBox::GridBox(Grid origin, Grid size) : _origin{origin}, _size{size} {}

Grid GridBox::origin() const { return _origin; }

Grid GridBox::size() const { return _size; }

Grid GridBox::center() const {
  return {_origin.i() + _size.i() / 2, _origin.j() + _size.j() / 2,
          _origin.k() + _size.k() / 2};
}

Grid GridBox::end() const { return origin() + size(); }

GridSpace::GridSpace(GridSpaceRegion region, double dx, double dy, double dz,
                     Dimension dimension, xt::xarray<double> e_node_x,
                     xt::xarray<double> e_node_y, xt::xarray<double> e_node_z)
    : _region{std::move(region)},
      _based_dx{dx},
      _based_dy{dy},
      _based_dz{dz},
      _dimension{dimension},
      _type{Type::UNIFORM},
      _e_node_x{std::move(e_node_x)},
      _e_node_y{std::move(e_node_y)},
      _e_node_z{std::move(e_node_z)} {}

GridSpace::Dimension GridSpace::dimension() const { return _dimension; }

GridSpace::Type GridSpace::type() const { return _type; }

GridSpace::GridSpaceRegion GridSpace::region() const { return _region; }

GridBox GridSpace::box() const {
  return GridBox{Grid{0, 0, 0}, Grid{sizeX(), sizeY(), sizeZ()}};
}

double GridSpace::basedDx() const { return _based_dx; }

double GridSpace::basedDy() const { return _based_dy; }

double GridSpace::basedDz() const { return _based_dz; }

double GridSpace::minDx() const { return _min_dx; }

double GridSpace::minDy() const { return _min_dy; }

double GridSpace::minDz() const { return _min_dz; }

const xt::xarray<double>& GridSpace::eNodeX() const { return _e_node_x; }

const xt::xarray<double>& GridSpace::eNodeY() const { return _e_node_y; }

const xt::xarray<double>& GridSpace::eNodeZ() const { return _e_node_z; }

const xt::xarray<double>& GridSpace::hNodeX() const { return _h_node_x; }

const xt::xarray<double>& GridSpace::hNodeY() const { return _h_node_y; }

const xt::xarray<double>& GridSpace::hNodeZ() const { return _h_node_z; }

const xt::xarray<double>& GridSpace::eSizeX() const { return _e_size_x; }

const xt::xarray<double>& GridSpace::eSizeY() const { return _e_size_y; }

const xt::xarray<double>& GridSpace::eSizeZ() const { return _e_size_z; }

const xt::xarray<double>& GridSpace::hSizeX() const { return _h_size_x; }

const xt::xarray<double>& GridSpace::hSizeY() const { return _h_size_y; }

const xt::xarray<double>& GridSpace::hSizeZ() const { return _h_size_z; }

const xt::xarray<double>& GridSpace::eSize(Axis::XYZ xyz) const {
  if (xyz == Axis::XYZ::X) {
    return eSizeX();
  }

  if (xyz == Axis::XYZ::Y) {
    return eSizeY();
  }

  if (xyz == Axis::XYZ::Z) {
    return eSizeZ();
  }

  throw XFDTDGridSpaceException{"Invalid xyz"};
}

const xt::xarray<double>& GridSpace::hSize(Axis::XYZ xyz) const {
  if (xyz == Axis::XYZ::X) {
    return hSizeX();
  }

  if (xyz == Axis::XYZ::Y) {
    return hSizeY();
  }

  if (xyz == Axis::XYZ::Z) {
    return hSizeZ();
  }

  throw XFDTDGridSpaceException{"Invalid xyz"};
}

std::size_t GridSpace::sizeX() const { return hNodeX().size(); }

std::size_t GridSpace::sizeY() const { return hNodeY().size(); }

std::size_t GridSpace::sizeZ() const { return hNodeZ().size(); }

Grid GridSpace::getGrid(double x, double y, double z) const {
  return {handleTransformX(x), handleTransformY(y), handleTransformZ(z)};
}

Grid GridSpace::getGrid(const Vector& vector) const {
  return getGrid(vector.x(), vector.y(), vector.z());
}

Vector GridSpace::getGridOriginVector(const Grid& grid) const {
  return Vector{_e_node_x(grid.i()), _e_node_y(grid.j()), _e_node_z(grid.k())};
}

Vector GridSpace::getGridCenterVector(const Grid& grid) const {
  return Vector{_h_node_x(grid.i()), _h_node_y(grid.j()), _h_node_z(grid.k())};
}

Vector GridSpace::getGridEndVector(const Grid& grid) const {
  return Vector{_e_node_x(grid.i() + 1), _e_node_y(grid.j() + 1),
                _e_node_z(grid.k() + 1)};
}

GridBox GridSpace::getGridBox(const Shape* shape) const {
  auto cube{shape->wrappedCube()};
  auto origin{getGrid(cube.origin())};
  auto end{getGrid(cube.end())};
  return GridBox{origin, end - origin};
}

void GridSpace::extendGridSpace(Axis::Direction direction, std::size_t num,
                                double dl) {
  switch (direction) {
    case Axis::Direction::XN: {  // insert e_node_x
      // auto front_offset{xt::arange<double>(eNodeX().front() - dl * num,
      //  eNodeX().front(), dl)};
      auto front_offset{eNodeX().front() - xt::arange<int>(num, 0, -1) * dl};
      eNodeX() = xt::concatenate(xt::xtuple(front_offset, eNodeX()), 0);
      break;
    }
    case Axis::Direction::XP: {  // insert e_node_x
      // auto back_offset{xt::arange<double>(
      // eNodeX().back() + dl, eNodeX().back() + dl * (num + 1), dl)};
      auto back_offset{eNodeX().back() + xt::arange<int>(1, num + 1) * dl};
      eNodeX() = xt::concatenate(xt::xtuple(eNodeX(), back_offset), 0);
      break;
    }
    case Axis::Direction::YN: {  // insert e_node_y
      auto front_offset{eNodeY().front() - xt::arange<int>(num, 0, -1) * dl};
      // auto front_offset{xt::arange<double>(eNodeY().front() - dl * num,
      //  eNodeY().front(), dl)};
      eNodeY() = xt::concatenate(xt::xtuple(front_offset, eNodeY()), 0);
      break;
    }
    case Axis::Direction::YP: {  // insert e_node_y
      // auto back_offset{xt::arange<double>(
      // eNodeY().back() + dl, eNodeY().back() + dl * (num + 1), dl)};
      auto back_offset{eNodeY().back() + xt::arange<int>(1, num + 1) * dl};
      eNodeY() = xt::concatenate(xt::xtuple(eNodeY(), back_offset), 0);
      break;
    }
    case Axis::Direction::ZN: {  // insert e_node_z
      // auto front_offset{xt::arange<double>(eNodeZ().front() - dl * num,
      //  eNodeZ().front(), dl)};
      auto front_offset{eNodeZ().front() - xt::arange<int>(num, 0, -1) * dl};
      eNodeZ() = xt::concatenate(xt::xtuple(front_offset, eNodeZ()), 0);
      break;
    }
    case Axis::Direction::ZP: {  // insert e_node_z
      // auto back_offset{xt::arange<double>(
      // eNodeZ().back() + dl, eNodeZ().back() + dl * (num + 1), dl)};
      auto back_offset{eNodeZ().back() + xt::arange<int>(1, num + 1) * dl};
      eNodeZ() = xt::concatenate(xt::xtuple(eNodeZ(), back_offset), 0);
      break;
    }
    default:
      throw XFDTDGridSpaceException{"Invalid direction"};
  };
}

xt::xarray<double>& GridSpace::eNodeX() { return _e_node_x; }

xt::xarray<double>& GridSpace::eNodeY() { return _e_node_y; }

xt::xarray<double>& GridSpace::eNodeZ() { return _e_node_z; }

xt::xarray<double>& GridSpace::hNodeX() { return _h_node_x; }

xt::xarray<double>& GridSpace::hNodeY() { return _h_node_y; }

xt::xarray<double>& GridSpace::hNodeZ() { return _h_node_z; }

xt::xarray<double>& GridSpace::eSizeX() { return _e_size_x; }

xt::xarray<double>& GridSpace::eSizeY() { return _e_size_y; }

xt::xarray<double>& GridSpace::eSizeZ() { return _e_size_z; }

xt::xarray<double>& GridSpace::hSizeX() { return _h_size_x; }

xt::xarray<double>& GridSpace::hSizeY() { return _h_size_y; }

xt::xarray<double>& GridSpace::hSizeZ() { return _h_size_z; }

xt::xarray<double>& GridSpace::eSize(Axis::XYZ xyz) {
  if (xyz == Axis::XYZ::X) {
    return eSizeX();
  }

  if (xyz == Axis::XYZ::Y) {
    return eSizeY();
  }

  if (xyz == Axis::XYZ::Z) {
    return eSizeZ();
  }

  throw XFDTDGridSpaceException{"Invalid xyz"};
}

xt::xarray<double>& GridSpace::hSize(Axis::XYZ xyz) {
  if (xyz == Axis::XYZ::X) {
    return hSizeX();
  }

  if (xyz == Axis::XYZ::Y) {
    return hSizeY();
  }

  if (xyz == Axis::XYZ::Z) {
    return hSizeZ();
  }

  throw XFDTDGridSpaceException{"Invalid xyz"};
}

void GridSpace::setMinDx(double min_dx) { _min_dx = min_dx; }

void GridSpace::setMinDy(double min_dy) { _min_dy = min_dy; }

void GridSpace::setMinDz(double min_dz) { _min_dz = min_dz; }

void GridSpace::generateGrid(std::size_t nx, std::size_t ny, std::size_t nz) {
  _grid.resize({nx, ny, nz});
  for (auto i{0}; i < nx; ++i) {
    for (auto j{0}; j < ny; ++j) {
      for (auto k{0}; k < nz; ++k) {
        _grid(i, j, k) = std::make_shared<Grid>(i, j, k);
      }
    }
  }
}

std::size_t GridSpace::handleTransformX(double x) const {
  if (x < _e_node_x.front()) {
    throw XFDTDGridSpaceException{"x is smaller than e_node_x.front()"};
  }

  if (x > _e_node_x.back()) {
    throw XFDTDGridSpaceException{"x is bigger than e_node_x.back()" +
                                  std::to_string(x) + " " +
                                  std::to_string(_e_node_x.back())};
  }

  return xt::argmin(xt::abs(_e_node_x - x)).front();
}

std::size_t GridSpace::handleTransformY(double y) const {
  if (y < _e_node_y.front()) {
    throw XFDTDGridSpaceException{"y is smaller than e_node_y.front()"};
  }

  if (y > _e_node_y.back() + basedDy() / 2) {
    throw XFDTDGridSpaceException{"y is bigger than e_node_y.back()"};
  }

  return xt::argmin(xt::abs(_e_node_y - y)).front();
}

std::size_t GridSpace::handleTransformZ(double z) const {
  if (z < _e_node_z.front()) {
    throw XFDTDGridSpaceException{"z is smaller than e_node_z.front()"};
  }

  if (!floatCompare(z, _e_node_z.back() + basedDz() / 2,
                    FloatCompareOperator::LessEqual)) {
    throw XFDTDGridSpaceException{"z is bigger than e_node_z.back() " +
                                  std::to_string(z) + " " +
                                  std::to_string(_e_node_z.back())};
  }

  return xt::argmin(xt::abs(_e_node_z - z)).front();
}

void GridSpace::correctGridSpaceForOne(double dl,
                                       const xt::xarray<double>& e_node,
                                       xt::xarray<double>& h_node,
                                       xt::xarray<double>& e_size,
                                       xt::xarray<double>& h_size) {
  h_node = calculateHNode(e_node);
  e_size = calculateESize(e_node);
  h_size = calculateHSize(h_node, e_node, dl);
}

xt::xarray<double> GridSpace::calculateHNode(const xt::xarray<double>& e_node) {
  return (xt::view(e_node, xt::range(_, -1)) +
          xt::view(e_node, xt::range(1, _))) /
         2;
}

xt::xarray<double> GridSpace::calculateESize(const xt::xarray<double>& e_node) {
  return (xt::view(e_node, xt::range(1, _)) -
          xt::view(e_node, xt::range(_, -1)));
}

xt::xarray<double> GridSpace::calculateHSize(const xt::xarray<double>& h_node,
                                             const xt::xarray<double>& e_node,
                                             double dl) {
  auto front_point{xt::xarray<double>{e_node.front() - dl / 2}};
  auto back_point{xt::xarray<double>{e_node.back() + dl / 2}};
  xt::xarray<double> temp{
      xt::concatenate(xt::xtuple(front_point, h_node, back_point))};
  return (xt::view(temp, xt::range(1, _)) - xt::view(temp, xt::range(_, -1)));
}

}  // namespace xfdtd
