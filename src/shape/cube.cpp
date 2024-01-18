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
  return std::string{"Cube("} + _origin.toString() + ", " + _size.toString() +
         ")";
}

Vector Cube::origin() const { return _origin; }

Vector Cube::size() const { return _size; }

Vector Cube::center() const { return _center; }

Vector Cube::end() const { return _end; }

double Cube::originX() const { return _origin.x(); }

double Cube::originY() const { return _origin.y(); }

double Cube::originZ() const { return _origin.z(); }

double Cube::sizeX() const { return _size.x(); }

double Cube::sizeY() const { return _size.y(); }

double Cube::sizeZ() const { return _size.z(); }

double Cube::centerX() const { return _center.x(); }

double Cube::centerY() const { return _center.y(); }

double Cube::centerZ() const { return _center.z(); }

double Cube::endX() const { return _end.x(); }

double Cube::endY() const { return _end.y(); }

double Cube::endZ() const { return _end.z(); }

bool Cube::isInside(double x, double y, double z) const {
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

Cube Cube::wrappedCube() const { return *this; }

void Cube::updateEnd() {
  _end = _origin + _size;
  if (std::isnan(_end.x())) {
    _end.setX(std::numeric_limits<double>::infinity());
  }
  if (std::isnan(_end.y())) {
    _end.setY(std::numeric_limits<double>::infinity());
  }
  if (std::isnan(_end.z())) {
    _end.setZ(std::numeric_limits<double>::infinity());
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
