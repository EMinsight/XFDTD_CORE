#include <xfdtd/shape/cube.h>
#include <xfdtd/shape/cylinder.h>
#include <xfdtd/shape/shape.h>

#include <utility>

#include "util/float_compare.h"
#include "xfdtd/coordinate_system/coordinate_system.h"

namespace xfdtd {

Cylinder::Cylinder(Vector center, double radius, double height, Axis axis)
    : _center{std::move(center)},
      _radius{radius},
      _height{height},
      _axis{std::move(axis)} {
  if (_radius < 0) {
    throw XFDTDShapeInvalidArgumentException{
        "Cylinder radius must be positive"};
  }
  if (_height < 0) {
    throw XFDTDShapeInvalidArgumentException{
        "Cylinder height must be positive"};
  }
}

std::unique_ptr<Shape> Cylinder::clone() const {
  return std::make_unique<Cylinder>(*this);
}

std::string Cylinder::toString() const {
  return std::string{"Cylinder("} + _center.toString() + ", " +
         std::to_string(_radius) + ", " + std::to_string(_height) + ", " +
         _axis.toString() + ")";
}

Vector Cylinder::center() const { return _center; }

double Cylinder::radius() const { return _radius; }

double Cylinder::height() const { return _height; }

Axis Cylinder::axis() const { return _axis; }

bool Cylinder::isInside(double x, double y, double z) const {
  auto wrapped_cube{wrappedCube()};
  if (!wrapped_cube.isInside(x, y, z)) {
    return false;
  }

  double dis{0};
  if (_axis == Axis::XYZ::X) {
    dis =
        std::sqrt(std::pow(y - _center.y(), 2) + std::pow(z - _center.z(), 2));
  } else if (_axis == Axis::XYZ::Y) {
    dis =
        std::sqrt(std::pow(x - _center.x(), 2) + std::pow(z - _center.z(), 2));
  } else if (_axis == Axis::XYZ::Z) {
    dis =
        std::sqrt(std::pow(x - _center.x(), 2) + std::pow(y - _center.y(), 2));
  } else {
    throw XFDTDShapeFailToSupportException{"Cylinder axis is not supported"};
  }

  return floatCompare(dis, radius(), FloatCompareOperator::LessEqual);
}

bool Cylinder::isInside(const Vector& vector) const {
  return isInside(vector.x(), vector.y(), vector.z());
}

Cube Cylinder::wrappedCube() const {
  auto size{radius() * 2};

  if (_axis == Axis::XYZ::X) {
    auto x{center().x() - _height / 2};
    auto y{center().y() - radius()};
    auto z{center().z() - radius()};
    return Cube{Vector{x, y, z}, Vector{height(), size, size}};
  }
  if (_axis == Axis::XYZ::Y) {
    auto x{center().x() - radius()};
    auto y{center().y() - _height / 2};
    auto z{center().z() - radius()};
    return Cube{Vector{x, y, z}, Vector{size, height(), size}};
  }
  if (_axis == Axis::XYZ::Z) {
    auto x{center().x() - radius()};
    auto y{center().y() - radius()};
    auto z{center().z() - _height / 2};
    return Cube{Vector{x, y, z}, Vector{size, size, height()}};
  }

  throw XFDTDShapeFailToSupportException{"Cylinder axis is not supported"};
}

}  // namespace xfdtd
