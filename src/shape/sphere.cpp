#include <xfdtd/shape/sphere.h>

#include "util/float_compare.h"
#include "xfdtd/coordinate_system/coordinate_system.h"
#include "xfdtd/shape/cube.h"

namespace xfdtd {

Sphere::Sphere(Vector center, Real radius)
    : _center{std::move(center)}, _radius{radius} {
  if (_radius < 0) {
    throw XFDTDShapeInvalidArgumentException{"Sphere radius must be positive"};
  }
}

std::unique_ptr<Shape> Sphere::clone() const {
  return std::make_unique<Sphere>(*this);
}

std::string Sphere::toString() const {
  return std::string{"Sphere("} + _center.toString() + ", " +
         std::to_string(_radius) + ")";
}

bool Sphere::isInside(Real x, Real y, Real z) const {
  return isInside(Vector{x, y, z});
}

bool Sphere::isInside(const Vector& vector) const {
  return floatCompare((vector - center()).normL2(), radius(),
                      FloatCompareOperator::LessEqual);
}

std::unique_ptr<Cube> Sphere::wrappedCube() const {
  return std::make_unique<Cube>(_center - Vector{_radius, _radius, _radius},
                                Vector{_radius * 2, _radius * 2, _radius * 2});
}

Vector Sphere::center() const { return _center; }

Real Sphere::radius() const { return _radius; }

Real Sphere::centerX() const { return _center.x(); }

Real Sphere::centerY() const { return _center.y(); }

Real Sphere::centerZ() const { return _center.z(); }

}  // namespace xfdtd
