#ifndef _XFDTD_LIB_SHAPE_H_
#define _XFDTD_LIB_SHAPE_H_

#include <xfdtd/coordinate_system/coordinate_system.h>

#include <memory>

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

  virtual bool isInside(double x, double y, double z) const = 0;

  virtual bool isInside(const Vector& vector) const = 0;

  virtual Cube wrappedCube() const = 0;
};

}  // namespace xfdtd

#endif  // _XFDTD_LIB_SHAPE_H_
