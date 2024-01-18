#include <xfdtd/shape/sphere.h>

#include "util/float_compare.h"
#include "xfdtd/coordinate_system/coordinate_system.h"
#include "xfdtd/shape/cube.h"

namespace xfdtd {

Sphere::Sphere(Vector center, double radius)
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

bool Sphere::isInside(double x, double y, double z) const {
  return isInside(Vector{x, y, z});
}

bool Sphere::isInside(const Vector& vector) const {
  return floatCompare((vector - center()).normL2(), radius(),
                      FloatCompareOperator::LessEqual);
}

Cube Sphere::wrappedCube() const {
  return Cube{_center - Vector{_radius, _radius, _radius},
              Vector{_radius * 2, _radius * 2, _radius * 2}};
}

Vector Sphere::center() const { return _center; }

double Sphere::radius() const { return _radius; }

double Sphere::centerX() const { return _center.x(); }

double Sphere::centerY() const { return _center.y(); }

double Sphere::centerZ() const { return _center.z(); }

}  // namespace xfdtd
