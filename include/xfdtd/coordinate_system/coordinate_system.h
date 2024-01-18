#ifndef _XFDTD_LIB_COORDINATE_SYSTEM_H_
#define _XFDTD_LIB_COORDINATE_SYSTEM_H_

#include <xfdtd/exception/exception.h>

#include <string>
#include <xtensor/xfixed.hpp>

namespace xfdtd {

class XFDTDCoordinateSystemException : public XFDTDException {
 public:
  explicit XFDTDCoordinateSystemException(
      std::string message = "XFDTD CoordinateSystem Exception")
      : XFDTDException{std::move(message)} {}
};

class XFDTDCoordinateSystemVectorException
    : public XFDTDCoordinateSystemException {
 public:
  explicit XFDTDCoordinateSystemVectorException(
      std::string message = "XFDTD CoordinateSystem Vector Exception")
      : XFDTDCoordinateSystemException{std::move(message)} {}
};

class XFDTDCoordinateSystemVectorNormalizedException
    : public XFDTDCoordinateSystemVectorException {
 public:
  explicit XFDTDCoordinateSystemVectorNormalizedException(
      std::string message = "XFDTD CoordinateSystem Vector Normalize Exception")
      : XFDTDCoordinateSystemVectorException{std::move(message)} {}
};

class XFDTDCoordinateSystemAxisException
    : public XFDTDCoordinateSystemException {
 public:
  explicit XFDTDCoordinateSystemAxisException(
      std::string message = "XFDTD CoordinateSystem Axis Exception")
      : XFDTDCoordinateSystemException{std::move(message)} {}
};

class XFDTDCoordinateSystemAxisDirectionException
    : public XFDTDCoordinateSystemAxisException {
 public:
  explicit XFDTDCoordinateSystemAxisDirectionException(
      std::string message = "XFDTD CoordinateSystem Axis Direction Exception")
      : XFDTDCoordinateSystemAxisException{std::move(message)} {}
};

using VectorData = xt::xtensor_fixed<double, xt::xshape<3>>;

class Vector {
 public:
  Vector() = default;

  explicit Vector(VectorData vector);

  Vector(double x, double y, double z);

  bool operator==(const Vector& rhs) const;

  bool operator!=(const Vector& rhs) const;

  Vector operator+(const Vector& rhs) const;

  Vector operator-(const Vector& rhs) const;

  Vector operator*(double rhs) const;

  Vector operator/(double rhs) const;

  virtual std::string toString() const;

  static Vector fromSpherical(double r, double theta, double phi);

  static Vector fromCylindrical(double r, double theta, double z);

  static Vector normalized(const Vector& vector);

  const VectorData& data() const;

  VectorData& data();

  double x() const;

  double y() const;

  double z() const;

  double normL2() const;

  void set(double x, double y, double z);

  void setX(double x);

  void setY(double y);

  void setZ(double z);

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

  static XYZ formDirectionToXYZ(Direction direction);

  static bool directionNegative(Axis::Direction direction);

  static bool directionPositive(Axis::Direction direction);
};

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

#endif  // _XFDTD_LIB_COORDINATE_SYSTEM_H_
