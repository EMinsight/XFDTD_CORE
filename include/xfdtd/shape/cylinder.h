#ifndef _XFDTD_CORE_CYLINDER_H_
#define _XFDTD_CORE_CYLINDER_H_

#include <xfdtd/coordinate_system/coordinate_system.h>
#include <xfdtd/shape/shape.h>

namespace xfdtd {

class Cylinder : public Shape {
 public:
  Cylinder(Vector center, Real radius, Real height, Axis axis);

  Cylinder(const Cylinder &) = default;

  Cylinder(Cylinder &&) noexcept = default;

  Cylinder &operator=(const Cylinder &) = default;

  Cylinder &operator=(Cylinder &&) noexcept = default;

  ~Cylinder() override = default;

  std::unique_ptr<Shape> clone() const override;

  std::string toString() const override;

  Vector center() const;

  Real radius() const;

  Real height() const;

  Axis axis() const;

  bool isInside(Real x, Real y, Real z) const override;

  bool isInside(const Vector& vector) const override;

  std::unique_ptr<Cube> wrappedCube() const override;

 private:
  Vector _center;
  Real _radius;
  Real _height;
  Axis _axis;
};

}  // namespace xfdtd

#endif  // _XFDTD_CORE_CYLINDER_H_
