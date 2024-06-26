#ifndef _XFDTD_CORE_COORDINATE_SYSTEM_H_
#define _XFDTD_CORE_COORDINATE_SYSTEM_H_

#include <xfdtd/common/type_define.h>
#include <xfdtd/exception/exception.h>

#include <string>
#include <xtensor/xfixed.hpp>

namespace xfdtd {

class XFDTDCoordinateSystemException : public XFDTDException {
 public:
  explicit XFDTDCoordinateSystemException(std::string message)
      : XFDTDException{std::move(message)} {}
};

class XFDTDCoordinateSystemVectorException
    : public XFDTDCoordinateSystemException {
 public:
  explicit XFDTDCoordinateSystemVectorException(std::string message)
      : XFDTDCoordinateSystemException{std::move(message)} {}
};

class XFDTDCoordinateSystemAxisException
    : public XFDTDCoordinateSystemException {
 public:
  explicit XFDTDCoordinateSystemAxisException(std::string message)
      : XFDTDCoordinateSystemException{std::move(message)} {}
};

class Vector {
 public:
  using VectorData = xt::xtensor_fixed<Real, xt::xshape<3>>;

 public:
  Vector() = default;

  explicit Vector(VectorData vector);

  Vector(Real x, Real y, Real z);

  bool operator==(const Vector& rhs) const;

  bool operator!=(const Vector& rhs) const;

  Vector operator+(const Vector& rhs) const;

  Vector operator-(const Vector& rhs) const;

  Vector operator*(Real rhs) const;

  Vector operator/(Real rhs) const;

  virtual std::string toString() const;

  static Vector fromSpherical(Real r, Real theta, Real phi);

  static Vector fromCylindrical(Real r, Real theta, Real z);

  static Vector normalized(const Vector& vector);

  const VectorData& data() const;

  VectorData& data();

  Real x() const;

  Real y() const;

  Real z() const;

  Real normL2() const;

  auto dot(const Vector& rhs) const -> Real;

  auto cross(const Vector& rhs) const -> Vector;

  void set(Real x, Real y, Real z);

  void setX(Real x);

  void setY(Real y);

  void setZ(Real z);

 private:
  VectorData _data{0, 0, 0};
};

class Axis : public Vector {
 public:
  enum class Direction {
    XN,
    XP,
    YN,
    YP,
    ZN,
    ZP,
  };

  enum class XYZ {
    X,
    Y,
    Z,
  };

  explicit Axis(const Vector& direction);

  bool operator==(const Axis& rhs) const;

  bool operator!=(const Axis& rhs) const;

  bool operator==(const Vector& rhs) const;

  bool operator!=(const Vector& rhs) const;

  bool operator==(Axis::Direction rhs) const;

  bool operator!=(Axis::Direction rhs) const;

  bool operator==(Axis::XYZ rhs) const;

  bool operator!=(Axis::XYZ rhs) const;

  std::string toString() const override;

  const static Axis XN;
  const static Axis XP;
  const static Axis YN;
  const static Axis YP;
  const static Axis ZN;
  const static Axis ZP;

  static Axis fromDirectionToAxis(Direction direction);

  static auto fromDirectionToXYZ(Direction direction) -> Axis::XYZ;

  static constexpr auto tangentialAAxis(Axis::XYZ c) -> Axis::XYZ;

  static constexpr auto tangentialBAxis(Axis::XYZ c) -> Axis::XYZ;

  static auto crossProduct(Axis::XYZ a, Axis::XYZ b) -> Axis::Direction;

  static bool directionNegative(Axis::Direction direction);

  static bool directionPositive(Axis::Direction direction);

  static auto toString(Axis::Direction direction) -> std::string;

  static auto toString(Axis::XYZ xyz) -> std::string;

  template <Axis::Direction direction>
  static constexpr auto fromDirectionToXYZ() -> Axis::XYZ;

  template <Axis::XYZ c>
  static constexpr auto tangentialAAxis() -> Axis::XYZ;

  template <Axis::XYZ c>
  static constexpr auto tangentialBAxis() -> Axis::XYZ;

  template <Axis::XYZ a, Axis::XYZ b>
  static constexpr auto crossProduct() -> Axis::Direction;

  template <Axis::Direction direction>
  static constexpr auto directionNegative() -> bool;

  template <Axis::Direction direction>
  static constexpr auto directionPositive() -> bool;
};

template <Axis::Direction direction>
inline constexpr auto Axis::fromDirectionToXYZ() -> Axis::XYZ {
  if constexpr (direction == Axis::Direction::XN ||
                direction == Axis::Direction::XP) {
    return Axis::XYZ::X;
  } else if constexpr (direction == Axis::Direction::YN ||
                       direction == Axis::Direction::YP) {
    return Axis::XYZ::Y;
  } else if constexpr (direction == Axis::Direction::ZN ||
                       direction == Axis::Direction::ZP) {
    return Axis::XYZ::Z;
  } else {
    throw XFDTDCoordinateSystemAxisException{
        "fromDirectionToXYZ: Invalid Axis::Direction value"};
  }
}

template <Axis::XYZ c>
inline constexpr auto Axis::tangentialAAxis() -> Axis::XYZ {
  if constexpr (c == Axis::XYZ::X) {
    return Axis::XYZ::Y;
  } else if constexpr (c == Axis::XYZ::Y) {
    return Axis::XYZ::Z;
  } else if constexpr (c == Axis::XYZ::Z) {
    return Axis::XYZ::X;
  } else {
    throw XFDTDCoordinateSystemAxisException{
        "tangentialAAxis: Invalid Axis::XYZ value"};
  }
}

template <Axis::XYZ c>
inline constexpr auto Axis::tangentialBAxis() -> Axis::XYZ {
  if constexpr (c == Axis::XYZ::X) {
    return Axis::XYZ::Z;
  } else if constexpr (c == Axis::XYZ::Y) {
    return Axis::XYZ::X;
  } else if constexpr (c == Axis::XYZ::Z) {
    return Axis::XYZ::Y;
  } else {
    throw XFDTDCoordinateSystemAxisException{
        "tangentialBAxis: Invalid Axis::XYZ value"};
  }
}

template <Axis::XYZ a, Axis::XYZ b>
inline constexpr auto Axis::crossProduct() -> Axis::Direction {
  if constexpr (a == Axis::XYZ::X && b == Axis::XYZ::Y) {
    return Axis::Direction::ZP;
  } else if constexpr (a == Axis::XYZ::Y && b == Axis::XYZ::Z) {
    return Axis::Direction::XP;
  } else if constexpr (a == Axis::XYZ::Z && b == Axis::XYZ::X) {
    return Axis::Direction::YP;
  } else if constexpr (a == Axis::XYZ::Y && b == Axis::XYZ::X) {
    return Axis::Direction::ZN;
  } else if constexpr (a == Axis::XYZ::Z && b == Axis::XYZ::Y) {
    return Axis::Direction::XN;
  } else if constexpr (a == Axis::XYZ::X && b == Axis::XYZ::Z) {
    return Axis::Direction::YN;
  } else {
    throw XFDTDCoordinateSystemAxisException{
        "crossProduct: Invalid Axis::XYZ value"};
  }
}

template <Axis::Direction direction>
inline constexpr auto Axis::directionNegative() -> bool {
  if constexpr (direction == Axis::Direction::XN ||
                direction == Axis::Direction::YN ||
                direction == Axis::Direction::ZN) {
    return true;
  } else if constexpr (direction == Axis::Direction::XP ||
                       direction == Axis::Direction::YP ||
                       direction == Axis::Direction::ZP) {
    return false;
  } else {
    throw XFDTDCoordinateSystemAxisException{
        "directionNegative: Invalid direction value"};
  }
}

template <Axis::Direction direction>
inline constexpr auto Axis::directionPositive() -> bool {
  if constexpr (direction == Axis::Direction::XN ||
                direction == Axis::Direction::YN ||
                direction == Axis::Direction::ZN) {
    return false;
  } else if constexpr (direction == Axis::Direction::XP ||
                       direction == Axis::Direction::YP ||
                       direction == Axis::Direction::ZP) {
    return true;
  } else {
    throw XFDTDCoordinateSystemAxisException{
        "directionPositive: Invalid direction value"};
  }
}

class CoordinateSystem {
 public:
  CoordinateSystem() = default;

  CoordinateSystem(Vector origin, Axis x, Axis y, Axis z);

  Vector origin() const;

  Axis x() const;

  Axis y() const;

  Axis z() const;

 private:
  Vector _origin{0, 0, 0};
  Axis _x{Axis::XP};
  Axis _y{Axis::YP};
  Axis _z{Axis::ZP};
};

}  // namespace xfdtd

#endif  // _XFDTD_CORE_COORDINATE_SYSTEM_H_
