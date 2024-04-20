#ifndef _XFDTD_CORE_SPHERE_H_
#define _XFDTD_CORE_SPHERE_H_

#include <xfdtd/shape/shape.h>

namespace xfdtd {

class Sphere : public Shape {
 public:
  Sphere(Vector center, Real radius);

  Sphere(const Sphere &) = default;

  Sphere(Sphere &&) noexcept = default;

  Sphere &operator=(const Sphere &) = default;

  Sphere &operator=(Sphere &&) noexcept = default;

  ~Sphere() override = default;

  std::unique_ptr<Shape> clone() const override;

  std::string toString() const override;

  bool isInside(Real x, Real y, Real z) const override;

  bool isInside(const Vector &vector) const override;

  std::unique_ptr<Cube> wrappedCube() const override;

  Vector center() const;

  Real radius() const;

  Real centerX() const;

  Real centerY() const;

  Real centerZ() const;

 private:
  Vector _center;
  Real _radius;
};

}  // namespace xfdtd

#endif  // _XFDTD_CORE_SPHERE_H_
