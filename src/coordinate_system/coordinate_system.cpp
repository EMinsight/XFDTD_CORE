#include <xfdtd/coordinate_system/coordinate_system.h>

#include <xtensor-blas/xlinalg.hpp>
#include <xtensor/xnorm.hpp>

namespace xfdtd {

CoordinateSystem::CoordinateSystem(Vector origin, Axis x, Axis y, Axis z)
    : _origin{std::move(origin)},
      _x{std::move(x)},
      _y{std::move(y)},
      _z{std::move(z)} {}

Vector CoordinateSystem::origin() const { return _origin; }

Axis CoordinateSystem::x() const { return _x; }

Axis CoordinateSystem::y() const { return _y; }

Axis CoordinateSystem::z() const { return _z; }

}  // namespace xfdtd
