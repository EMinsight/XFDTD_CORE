#ifndef _XFDTD_CORE_SHAPE_H_
#define _XFDTD_CORE_SHAPE_H_

#include <xfdtd/common/type_define.h>
#include <xfdtd/coordinate_system/coordinate_system.h>

#include <memory>
#include <vector>

namespace xfdtd {

class XFDTDShapeException : public XFDTDException {
 public:
  XFDTDShapeException() = default;

  explicit XFDTDShapeException(const std::string& message)
      : XFDTDException{message} {}
};

class XFDTDShapeInvalidArgumentException : public XFDTDShapeException {
 public:
  XFDTDShapeInvalidArgumentException() = default;

  explicit XFDTDShapeInvalidArgumentException(const std::string& message)
      : XFDTDShapeException{message} {}
};

class XFDTDShapeFailToSupportException : public XFDTDShapeException {
 public:
  XFDTDShapeFailToSupportException() = default;

  explicit XFDTDShapeFailToSupportException(const std::string& message)
      : XFDTDShapeException{message} {}
};

class Cube;

class Shape {
 public:
  Shape() = default;

  Shape(const Shape&) = default;

  Shape(Shape&&) noexcept = default;

  Shape& operator=(const Shape&) = default;

  Shape& operator=(Shape&&) noexcept = default;

  virtual ~Shape() = default;

  virtual std::unique_ptr<Shape> clone() const = 0;

  virtual std::string toString() const;

  virtual bool isInside(Real x, Real y, Real z, Real eps) const = 0;

  virtual bool isInside(const Vector& vector, Real eps) const = 0;

  virtual std::unique_ptr<Cube> wrappedCube() const = 0;

  static std::unique_ptr<Cube> makeWrappedCube(
      const std::vector<const Shape*>& shapes);
};

}  // namespace xfdtd

#endif  // _XFDTD_CORE_SHAPE_H_
