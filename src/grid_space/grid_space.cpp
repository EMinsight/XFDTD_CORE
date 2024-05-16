#include <xfdtd/coordinate_system/coordinate_system.h>
#include <xfdtd/grid_space/grid_space.h>

#include <algorithm>
#include <memory>
#include <sstream>
#include <utility>
#include <xtensor.hpp>

#include "util/float_compare.h"
#include "xfdtd/common/type_define.h"
#include "xfdtd/shape/cube.h"
#include "xfdtd/shape/shape.h"

namespace xfdtd {

Grid::Grid(std::size_t i, std::size_t j, std::size_t k,
           std::size_t material_index)
    : _i{i}, _j{j}, _k{k}, _material_index{material_index} {}

Grid Grid::operator+(const Grid& grid) const {
  return {_i + grid._i, _j + grid._j, _k + grid._k};
}

Grid Grid::operator-(const Grid& grid) const {
  return {_i - grid._i, _j - grid._j, _k - grid._k};
}

std::size_t Grid::i() const { return _i; }

std::size_t Grid::j() const { return _j; }

std::size_t Grid::k() const { return _k; }

std::size_t Grid::materialIndex() const { return _material_index; }

void Grid::setMaterialIndex(std::size_t index) { _material_index = index; }

std::string Grid::toString() const {
  return "Grid{i: " + std::to_string(_i) + ", j: " + std::to_string(_j) +
         ", k: " + std::to_string(_k) + "}";
}

GridBox::GridBox(Grid origin, Grid size) : _origin{origin}, _size{size} {}

Grid GridBox::origin() const { return _origin; }

Grid GridBox::size() const { return _size; }

Grid GridBox::center() const {
  return {_origin.i() + _size.i() / 2, _origin.j() + _size.j() / 2,
          _origin.k() + _size.k() / 2};
}

Grid GridBox::end() const { return origin() + size(); }

std::string GridBox::toString() const {
  std::stringstream ss;
  ss << "GridBox:\n";
  ss << " Origin: " << _origin.toString() << "\n";
  ss << " Size: " << _size.toString();
  return ss.str();
}

GridSpace::GridSpace(Real dx, Real dy, Real dz, Dimension dimension,
                     Array1D<Real> e_node_x, Array1D<Real> e_node_y,
                     Array1D<Real> e_node_z)
    : _based_dx{dx},
      _based_dy{dy},
      _based_dz{dz},
      _dimension{dimension},
      _type{Type::UNIFORM},
      _e_node_x{std::move(e_node_x)},
      _e_node_y{std::move(e_node_y)},
      _e_node_z{std::move(e_node_z)} {}

GridSpace::GridSpace(Dimension dimension, Type type, GridBox global_box,
                     Real based_dx, Real based_dy, Real based_dz, Real min_dx,
                     Real min_dy, Real min_dz, Array1D<Real> e_node_x,
                     Array1D<Real> e_node_y, Array1D<Real> e_node_z,
                     Array1D<Real> h_node_x, Array1D<Real> h_node_y,
                     Array1D<Real> h_node_z, Array1D<Real> e_size_x,
                     Array1D<Real> e_size_y, Array1D<Real> e_size_z,
                     Array1D<Real> h_size_x, Array1D<Real> h_size_y,
                     Array1D<Real> h_size_z)
    : _dimension{dimension},
      _type{type},
      _global_box{global_box},
      _based_dx{based_dx},
      _based_dy{based_dy},
      _based_dz{based_dz},
      _min_dx{min_dx},
      _min_dy{min_dy},
      _min_dz{min_dz},
      _e_node_x{std::move(e_node_x)},
      _e_node_y{std::move(e_node_y)},
      _e_node_z{std::move(e_node_z)},
      _h_node_x{std::move(h_node_x)},
      _h_node_y{std::move(h_node_y)},
      _h_node_z{std::move(h_node_z)},
      _e_size_x{std::move(e_size_x)},
      _e_size_y{std::move(e_size_y)},
      _e_size_z{std::move(e_size_z)},
      _h_size_x{std::move(h_size_x)},
      _h_size_y{std::move(h_size_y)},
      _h_size_z{std::move(h_size_z)} {}

GridSpace::Dimension GridSpace::dimension() const { return _dimension; }

GridSpace::Type GridSpace::type() const { return _type; }

Cube GridSpace::region() const {
  return Cube{Vector{eNodeX().front(), eNodeY().front(), eNodeZ().front()},
              Vector{eNodeX().back() - eNodeX().front(),
                     eNodeY().back() - eNodeY().front(),
                     eNodeZ().back() - eNodeZ().front()}};
}

GridBox GridSpace::box() const {
  return GridBox{Grid{0, 0, 0}, Grid{sizeX(), sizeY(), sizeZ()}};
}

Real GridSpace::basedDx() const { return _based_dx; }

Real GridSpace::basedDy() const { return _based_dy; }

Real GridSpace::basedDz() const { return _based_dz; }

Real GridSpace::minDx() const { return _min_dx; }

Real GridSpace::minDy() const { return _min_dy; }

Real GridSpace::minDz() const { return _min_dz; }

std::shared_ptr<GridSpace> GridSpace::globalGridSpace() const {
  return _global_grid_space.lock();
}

const Array1D<Real>& GridSpace::eNodeX() const { return _e_node_x; }

const Array1D<Real>& GridSpace::eNodeY() const { return _e_node_y; }

const Array1D<Real>& GridSpace::eNodeZ() const { return _e_node_z; }

const Array1D<Real>& GridSpace::hNodeX() const { return _h_node_x; }

const Array1D<Real>& GridSpace::hNodeY() const { return _h_node_y; }

const Array1D<Real>& GridSpace::hNodeZ() const { return _h_node_z; }

const Array1D<Real>& GridSpace::eSizeX() const { return _e_size_x; }

const Array1D<Real>& GridSpace::eSizeY() const { return _e_size_y; }

const Array1D<Real>& GridSpace::eSizeZ() const { return _e_size_z; }

const Array1D<Real>& GridSpace::hSizeX() const { return _h_size_x; }

const Array1D<Real>& GridSpace::hSizeY() const { return _h_size_y; }

const Array1D<Real>& GridSpace::hSizeZ() const { return _h_size_z; }

const Array1D<Real>& GridSpace::eSize(Axis::XYZ xyz) const {
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

const Array1D<Real>& GridSpace::hSize(Axis::XYZ xyz) const {
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

Grid GridSpace::getGrid(Real x, Real y, Real z) const {
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
  auto origin{getGrid(cube->origin())};
  auto end{getGrid(cube->end())};
  return GridBox{origin, end - origin};
}

void GridSpace::extendGridSpace(Axis::Direction direction, std::size_t num,
                                Real dl) {
  switch (direction) {
    case Axis::Direction::XN: {  // insert e_node_x
      // auto front_offset{xt::arange<Real>(eNodeX().front() - dl * num,
      //  eNodeX().front(), dl)};
      auto front_offset{eNodeX().front() - xt::arange<int>(num, 0, -1) * dl};
      eNodeX() = xt::concatenate(xt::xtuple(front_offset, eNodeX()), 0);
      break;
    }
    case Axis::Direction::XP: {  // insert e_node_x
      // auto back_offset{xt::arange<Real>(
      // eNodeX().back() + dl, eNodeX().back() + dl * (num + 1), dl)};
      auto back_offset{eNodeX().back() + xt::arange<int>(1, num + 1) * dl};
      eNodeX() = xt::concatenate(xt::xtuple(eNodeX(), back_offset), 0);
      break;
    }
    case Axis::Direction::YN: {  // insert e_node_y
      auto front_offset{eNodeY().front() - xt::arange<int>(num, 0, -1) * dl};
      // auto front_offset{xt::arange<Real>(eNodeY().front() - dl * num,
      //  eNodeY().front(), dl)};
      eNodeY() = xt::concatenate(xt::xtuple(front_offset, eNodeY()), 0);
      break;
    }
    case Axis::Direction::YP: {  // insert e_node_y
      // auto back_offset{xt::arange<Real>(
      // eNodeY().back() + dl, eNodeY().back() + dl * (num + 1), dl)};
      auto back_offset{eNodeY().back() + xt::arange<int>(1, num + 1) * dl};
      eNodeY() = xt::concatenate(xt::xtuple(eNodeY(), back_offset), 0);
      break;
    }
    case Axis::Direction::ZN: {  // insert e_node_z
      // auto front_offset{xt::arange<Real>(eNodeZ().front() - dl * num,
      //  eNodeZ().front(), dl)};
      auto front_offset{eNodeZ().front() - xt::arange<int>(num, 0, -1) * dl};
      eNodeZ() = xt::concatenate(xt::xtuple(front_offset, eNodeZ()), 0);
      break;
    }
    case Axis::Direction::ZP: {  // insert e_node_z
      // auto back_offset{xt::arange<Real>(
      // eNodeZ().back() + dl, eNodeZ().back() + dl * (num + 1), dl)};
      auto back_offset{eNodeZ().back() + xt::arange<int>(1, num + 1) * dl};
      eNodeZ() = xt::concatenate(xt::xtuple(eNodeZ(), back_offset), 0);
      break;
    }
    default:
      throw XFDTDGridSpaceException{"Invalid direction"};
  };
}

GridBox GridSpace::getGridBoxWithoutCheck(const Shape* shape) const {
  auto cube{shape->wrappedCube()};
  auto origin{getGridWithoutCheck(cube->origin())};
  auto end{getGridWithoutCheck(cube->end())};
  return GridBox{origin, end - origin};
}

Grid GridSpace::getGridWithoutCheck(const Vector& vector) const {
  return {handleTransformXWithoutCheck(vector.x()),
          handleTransformYWithoutCheck(vector.y()),
          handleTransformZWithoutCheck(vector.z())};
}

Array1D<Real>& GridSpace::eNodeX() { return _e_node_x; }

Array1D<Real>& GridSpace::eNodeY() { return _e_node_y; }

Array1D<Real>& GridSpace::eNodeZ() { return _e_node_z; }

Array1D<Real>& GridSpace::hNodeX() { return _h_node_x; }

Array1D<Real>& GridSpace::hNodeY() { return _h_node_y; }

Array1D<Real>& GridSpace::hNodeZ() { return _h_node_z; }

Array1D<Real>& GridSpace::eSizeX() { return _e_size_x; }

Array1D<Real>& GridSpace::eSizeY() { return _e_size_y; }

Array1D<Real>& GridSpace::eSizeZ() { return _e_size_z; }

Array1D<Real>& GridSpace::hSizeX() { return _h_size_x; }

Array1D<Real>& GridSpace::hSizeY() { return _h_size_y; }

Array1D<Real>& GridSpace::hSizeZ() { return _h_size_z; }

Array1D<Real>& GridSpace::eSize(Axis::XYZ xyz) {
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

Array1D<Real>& GridSpace::hSize(Axis::XYZ xyz) {
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

void GridSpace::setMinDx(Real min_dx) { _min_dx = min_dx; }

void GridSpace::setMinDy(Real min_dy) { _min_dy = min_dy; }

void GridSpace::setMinDz(Real min_dz) { _min_dz = min_dz; }

void GridSpace::generateMaterialGrid(std::size_t nx, std::size_t ny,
                                     std::size_t nz) {
  _grid_with_material =
      Array3D<std::shared_ptr<Grid>>::from_shape({nx, ny, nz});
  for (auto i{0}; i < nx; ++i) {
    for (auto j{0}; j < ny; ++j) {
      for (auto k{0}; k < nz; ++k) {
        _grid_with_material(i, j, k) = std::make_shared<Grid>(i, j, k);
      }
    }
  }
}

void GridSpace::setGlobalGridSpace(std::weak_ptr<GridSpace> global_grid_space) {
  auto temp{global_grid_space.lock()};
  if (temp == nullptr) {
    throw XFDTDGridSpaceException{"Invalid global_grid_space"};
  }

  if (temp->dimension() != dimension()) {
    throw XFDTDGridSpaceException{"Invalid dimension"};
  }

  _global_grid_space = std::move(global_grid_space);
}

std::size_t GridSpace::handleTransformX(Real x) const {
  if (x < _e_node_x.front()) {
    throw XFDTDGridSpaceException{"x is smaller than e_node_x.front()"};
  }

  if (!floatCompare(x, _e_node_x.back(), FloatCompareOperator::LessEqual, eps())) {
    throw XFDTDGridSpaceException{"x is bigger than e_node_x.back()" +
                                  std::to_string(x) + " " +
                                  std::to_string(_e_node_x.back())};
  }

  return xt::argmin(xt::abs(_e_node_x - x)).front();
}

GridBox GridSpace::validGridBoxEx() const {
  return {Grid{0, 1, 1}, Grid{sizeX(), sizeY() - 1, sizeZ() - 1}};
}

GridBox GridSpace::validGridBoxEy() const {
  return {Grid{1, 0, 1}, Grid{sizeX() - 1, sizeY(), sizeZ() - 1}};
}

GridBox GridSpace::validGridBoxEz() const {
  return {Grid{1, 1, 0}, Grid{sizeX() - 1, sizeY() - 1, sizeZ()}};
}

GridBox GridSpace::globalBox() const { return _global_box; }

static std::string dimensionToString(GridSpace::Dimension condition) {
  switch (condition) {
    case GridSpace::Dimension::ONE:
      return "ONE";
    case GridSpace::Dimension::TWO:
      return "TWO";
    case GridSpace::Dimension::THREE:
      return "THREE";
    default:
      return "UNDEFINED";
  }
}

static std::string typeToString(GridSpace::Type condition) {
  switch (condition) {
    case GridSpace::Type::UNIFORM:
      return "UNIFORM";
    case GridSpace::Type::NONUNIFORM:
      return "NON_UNIFORM";
    default:
      return "UNDEFINED";
  }
}

std::string GridSpace::toString() const {
  std::stringstream ss;

  ss << "GridSpace: \n";
  ss << " Dimension: " << dimensionToString(dimension()) << "\n";
  ss << " Type: " << typeToString(type()) << "\n";
  ss << " Region: " << region().toString() << "\n";
  ss << "  " << box().toString() << "\n";
  ss << " Global " << globalBox().toString();

  return ss.str();
}

Grid GridSpace::transformNodeToGlobal(const Grid& grid) const {
  auto offset{globalBox().origin()};
  return offset + grid;
}

GridBox GridSpace::transformNodeToGlobal(const GridBox& grid_box) const {
  auto offset{globalBox().origin()};
  return {offset + grid_box.origin(), grid_box.size()};
}

auto GridSpace::eps() const -> Real {
  constexpr auto num = 1000;
  if (dimension() == Dimension::ONE) {
    return minDz() / num;
  }
  if (dimension() == Dimension::TWO) {
    return std::min(minDx(), minDy()) / num;
  }
  return std::min({minDx(), minDy(), minDz()}) / num;
}

std::size_t GridSpace::handleTransformY(Real y) const {
  if (y < _e_node_y.front()) {
    throw XFDTDGridSpaceException{"y is smaller than e_node_y.front()"};
  }

  if (y > _e_node_y.back() + basedDy() / 2) {
    throw XFDTDGridSpaceException{"y is bigger than e_node_y.back()"};
  }

  return xt::argmin(xt::abs(_e_node_y - y)).front();
}

std::size_t GridSpace::handleTransformZ(Real z) const {
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

std::size_t GridSpace::handleTransformXWithoutCheck(Real x) const {
  return xt::argmin(xt::abs(_e_node_x - x)).front();
}

std::size_t GridSpace::handleTransformYWithoutCheck(Real y) const {
  return xt::argmin(xt::abs(_e_node_y - y)).front();
}

std::size_t GridSpace::handleTransformZWithoutCheck(Real z) const {
  return xt::argmin(xt::abs(_e_node_z - z)).front();
}

void GridSpace::correctGridSpaceForOne(Real dl, const Array1D<Real>& e_node,
                                       Array1D<Real>& h_node,
                                       Array1D<Real>& e_size,
                                       Array1D<Real>& h_size) {
  h_node = calculateHNode(e_node);
  e_size = calculateESize(e_node);
  h_size = calculateHSize(h_node, e_node, dl);
}

Array1D<Real> GridSpace::calculateHNode(const Array1D<Real>& e_node) {
  return (xt::view(e_node, xt::range(_, -1)) +
          xt::view(e_node, xt::range(1, _))) /
         2;
}

Array1D<Real> GridSpace::calculateESize(const Array1D<Real>& e_node) {
  return (xt::view(e_node, xt::range(1, _)) -
          xt::view(e_node, xt::range(_, -1)));
}

Array1D<Real> GridSpace::calculateHSize(const Array1D<Real>& h_node,
                                        const Array1D<Real>& e_node, Real dl) {
  auto front_point{Array1D<Real>{e_node.front() - dl / 2}};
  auto back_point{Array1D<Real>{e_node.back() + dl / 2}};
  Array1D<Real> temp{
      xt::concatenate(xt::xtuple(front_point, h_node, back_point))};
  return (xt::view(temp, xt::range(1, _)) - xt::view(temp, xt::range(_, -1)));
}

}  // namespace xfdtd
