#ifndef _XFDTD_LIB_CUBE_H_
#define _XFDTD_LIB_CUBE_H_

#include <xfdtd/coordinate_system/coordinate_system.h>
#include <xfdtd/shape/shape.h>

#include <memory>

namespace xfdtd {

class Cube : public Shape {
 public:
  Cube() = default;

  Cube(Vector origin, Vector size);

  Cube(const Cube&) = default;

  Cube(Cube&&) noexcept = default;

  Cube& operator=(const Cube&) = default;

  Cube& operator=(Cube&&) noexcept = default;

  ~Cube() override = default;

  std::unique_ptr<Shape> clone() const override;

  std::string toString() const override;

  Vector origin() const;

  Vector size() const;

  Vector center() const;

  Vector end() const;

  double originX() const;

  double originY() const;

  double originZ() const;

  double sizeX() const;

  double sizeY() const;

  double sizeZ() const;

  double centerX() const;

  double centerY() const;

  double centerZ() const;

  double endX() const;

  double endY() const;

  double endZ() const;

  bool isInside(double x, double y, double z) const override;

  bool isInside(const Vector& vector) const override;

  Cube wrappedCube() const override;

 private:
  Vector _origin{};
  Vector _size{};
  Vector _end{};
  Vector _center{};

  void updateEnd();
  void updateCenter();
};

}  // namespace xfdtd

#endif  // _XFDTD_LIB_CUBE_H_
