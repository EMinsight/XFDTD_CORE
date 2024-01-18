#include <xfdtd/coordinate_system/coordinate_system.h>

#include <cmath>
#include <string>
#include <utility>
#include <xtensor.hpp>

#include "util/float_compare.h"

namespace xfdtd {

Vector::Vector(xt::xtensor_fixed<double, xt::xshape<3>> vector)
    : _data{std::move(vector)} {}

Vector::Vector(double x, double y, double z) : _data{x, y, z} {}

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

Vector Vector::operator*(double rhs) const { return Vector{_data * rhs}; }

Vector Vector::operator/(double rhs) const { return Vector{_data / rhs}; }

std::string Vector::toString() const {
  return std::string{"Vector("} + std::to_string(x()) + ", " +
         std::to_string(y()) + ", " + std::to_string(z()) + ")";
}

Vector Vector::fromSpherical(double r, double theta, double phi) {
  return {r * std::sin(theta) * std::cos(phi),
          r * std::sin(theta) * std::sin(phi), r * std::cos(theta)};
}

Vector Vector::fromCylindrical(double r, double theta, double z) {
  return {r * std::cos(theta), r * std::sin(theta), z};
}

Vector Vector::normalized(const Vector& vector) {
  auto norm_l2{vector.normL2()};
  if (norm_l2 == 0) {
    // double equal to zero?
    throw XFDTDCoordinateSystemVectorNormalizedException{};
  }

  return Vector{vector.data() / norm_l2};
}

const VectorData& Vector::data() const { return _data; }

VectorData& Vector::data() { return _data; }

double Vector::x() const { return _data[0]; }

double Vector::y() const { return _data[1]; }

double Vector::z() const { return _data[2]; }

double Vector::normL2() const { return xt::eval(xt::norm_l2(_data))(); }

void Vector::set(double x, double y, double z) { _data = {x, y, z}; }

void Vector::setX(double x) { _data[0] = x; }

void Vector::setY(double y) { _data[1] = y; }

void Vector::setZ(double z) { _data[2] = z; }

}  // namespace xfdtd
