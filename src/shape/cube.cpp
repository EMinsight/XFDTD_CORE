#include <xfdtd/shape/cube.h>

#include <cmath>
#include <limits>
#include <utility>

#include "util/float_compare.h"

namespace xfdtd {

Cube::Cube(Vector origin, Vector size)
    : _origin{std::move(origin)}, _size{std::move(size)} {
  updateEnd();
  updateCenter();
}

std::unique_ptr<Shape> Cube::clone() const {
  return std::make_unique<Cube>(*this);
}

std::string Cube::toString() const {
  return std::string{"Cube(origin: "} + _origin.toString() +
         ", size:" + _size.toString() + " end: " + _end.toString() + ")";
}

Vector Cube::origin() const { return _origin; }

Vector Cube::size() const { return _size; }

Vector Cube::center() const { return _center; }

Vector Cube::end() const { return _end; }

Real Cube::originX() const { return _origin.x(); }

Real Cube::originY() const { return _origin.y(); }

Real Cube::originZ() const { return _origin.z(); }

Real Cube::sizeX() const { return _size.x(); }

Real Cube::sizeY() const { return _size.y(); }

Real Cube::sizeZ() const { return _size.z(); }

Real Cube::centerX() const { return _center.x(); }

Real Cube::centerY() const { return _center.y(); }

Real Cube::centerZ() const { return _center.z(); }

Real Cube::endX() const { return _end.x(); }

Real Cube::endY() const { return _end.y(); }

Real Cube::endZ() const { return _end.z(); }

bool Cube::isInside(Real x, Real y, Real z) const {
  auto x_inside{floatCompare(originX(), x, FloatCompareOperator::LessEqual) &&
                floatCompare(x, endX(), FloatCompareOperator::LessEqual)};
  auto y_inside{floatCompare(originY(), y, FloatCompareOperator::LessEqual) &&
                floatCompare(y, endY(), FloatCompareOperator::LessEqual)};
  auto z_inside{floatCompare(originZ(), z, FloatCompareOperator::LessEqual) &&
                floatCompare(z, endZ(), FloatCompareOperator::LessEqual)};
  return x_inside && y_inside && z_inside;
}

bool Cube::isInside(const Vector& vector) const {
  return isInside(vector.x(), vector.y(), vector.z());
}

std::unique_ptr<Cube> Cube::wrappedCube() const {
  return std::make_unique<Cube>(*this);
}

void Cube::updateEnd() {
  _end = _origin + _size;
  if (std::isnan(_end.x())) {
    _end.setX(std::numeric_limits<Real>::infinity());
  }
  if (std::isnan(_end.y())) {
    _end.setY(std::numeric_limits<Real>::infinity());
  }
  if (std::isnan(_end.z())) {
    _end.setZ(std::numeric_limits<Real>::infinity());
  }
}

void Cube::updateCenter() {
  _center = _origin + _size / 2;
  if (std::isnan(_center.x())) {
    _center.setX(0);
  }
  if (std::isnan(_center.y())) {
    _center.setY(0);
  }
  if (std::isnan(_center.z())) {
    _center.setZ(0);
  }
}

}  // namespace xfdtd
