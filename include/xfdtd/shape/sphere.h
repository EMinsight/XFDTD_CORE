#ifndef _XFDTD_CORE_SPHERE_H_
#define _XFDTD_CORE_SPHERE_H_

#include <xfdtd/shape/shape.h>

namespace xfdtd {

class Sphere : public Shape {
 public:
  Sphere(Vector center, double radius);

  Sphere(const Sphere &) = default;

  Sphere(Sphere &&) noexcept = default;

  Sphere &operator=(const Sphere &) = default;

  Sphere &operator=(Sphere &&) noexcept = default;

  ~Sphere() override = default;

  std::unique_ptr<Shape> clone() const override;

  std::string toString() const override;

  bool isInside(double x, double y, double z) const override;

  bool isInside(const Vector &vector) const override;

  std::unique_ptr<Cube> wrappedCube() const override;

  Vector center() const;

  double radius() const;

  double centerX() const;

  double centerY() const;

  double centerZ() const;

 private:
  Vector _center;
  double _radius;
};

}  // namespace xfdtd

#endif  // _XFDTD_CORE_SPHERE_H_
