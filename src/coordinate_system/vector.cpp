#include <xfdtd/coordinate_system/coordinate_system.h>

#include <cmath>
#include <string>
#include <utility>
#include <xtensor.hpp>

#include "util/float_compare.h"

namespace xfdtd {

Vector::Vector(xt::xtensor_fixed<Real, xt::xshape<3>> vector)
    : _data{std::move(vector)} {}

Vector::Vector(Real x, Real y, Real z) : _data{x, y, z} {}

bool Vector::operator==(const Vector& rhs) const {
  return floatCompareEqual(x(), rhs.x()) && floatCompareEqual(y(), rhs.y()) &&
         floatCompareEqual(z(), rhs.z());
}

bool Vector::operator!=(const Vector& rhs) const { return !(*this == rhs); }

Vector Vector::operator+(const Vector& rhs) const {
  return Vector{_data + rhs.data()};
}

Vector Vector::operator-(const Vector& rhs) const {
  return Vector{_data - rhs.data()};
}

Vector Vector::operator*(Real rhs) const { return Vector{_data * rhs}; }

Vector Vector::operator/(Real rhs) const { return Vector{_data / rhs}; }

std::string Vector::toString() const {
  return std::string{"Vector("} + std::to_string(x()) + ", " +
         std::to_string(y()) + ", " + std::to_string(z()) + ")";
}

Vector Vector::fromSpherical(Real r, Real theta, Real phi) {
  return {r * std::sin(theta) * std::cos(phi),
          r * std::sin(theta) * std::sin(phi), r * std::cos(theta)};
}

Vector Vector::fromCylindrical(Real r, Real theta, Real z) {
  return {r * std::cos(theta), r * std::sin(theta), z};
}

Vector Vector::normalized(const Vector& vector) {
  auto norm_l2{vector.normL2()};
  if (norm_l2 == 0) {
    // Real equal to zero?
    throw XFDTDCoordinateSystemVectorException{
        "Cannot normalize a vector with zero norm."};
  }

  return Vector{vector.data() / norm_l2};
}

const Vector::VectorData& Vector::data() const { return _data; }

Vector::VectorData& Vector::data() { return _data; }

Real Vector::x() const { return _data[0]; }

Real Vector::y() const { return _data[1]; }

Real Vector::z() const { return _data[2]; }

Real Vector::normL2() const { return xt::eval(xt::norm_l2(_data))(); }

auto Vector::dot(const Vector& rhs) const -> Real {
  return x() * rhs.x() + y() * rhs.y() + z() * rhs.z();
}

auto Vector::cross(const Vector& rhs) const -> Vector {
  return {y() * rhs.z() - z() * rhs.y(), z() * rhs.x() - x() * rhs.z(),
          x() * rhs.y() - y() * rhs.x()};
}

auto Vector::orthogonal(const Vector& rhs) const -> bool {
  return floatCompareEqual(dot(rhs), 0);
}

void Vector::set(Real x, Real y, Real z) { _data = {x, y, z}; }

void Vector::setX(Real x) { _data[0] = x; }

void Vector::setY(Real y) { _data[1] = y; }

void Vector::setZ(Real z) { _data[2] = z; }

}  // namespace xfdtd
