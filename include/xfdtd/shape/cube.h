#ifndef _XFDTD_CORE_CUBE_H_
#define _XFDTD_CORE_CUBE_H_

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

  Real originX() const;

  Real originY() const;

  Real originZ() const;

  Real sizeX() const;

  Real sizeY() const;

  Real sizeZ() const;

  Real centerX() const;

  Real centerY() const;

  Real centerZ() const;

  Real endX() const;

  Real endY() const;

  Real endZ() const;

  bool isInside(Real x, Real y, Real z) const override;

  bool isInside(const Vector& vector) const override;

  std::unique_ptr<Cube> wrappedCube() const override;

 private:
  Vector _origin{};
  Vector _size{};
  Vector _end{};
  Vector _center{};

  void updateEnd();
  void updateCenter();
};

}  // namespace xfdtd

#endif  // _XFDTD_CORE_CUBE_H_
