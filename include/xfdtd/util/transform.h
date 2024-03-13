#ifndef _XFDTD_CORE_TRANSFORM_H_
#define _XFDTD_CORE_TRANSFORM_H_

#include <tuple>

#include "xfdtd/coordinate_system/coordinate_system.h"

namespace xfdtd::transform {

template <typename T>
inline auto xYZToABC(const std::tuple<T, T, T> &data, Axis::XYZ xyz) {
  switch (xyz) {
    case Axis::XYZ::X: {
      // a: y b: z c: x
      return std::make_tuple(std::get<1>(data), std::get<2>(data),
                             std::get<0>(data));
    }
    case Axis::XYZ::Y: {
      // a: z b: x c: y
      return std::make_tuple(std::get<2>(data), std::get<0>(data),
                             std::get<1>(data));
    }
    case Axis::XYZ::Z: {
      // a: x b: y c: z
      return std::make_tuple(std::get<0>(data), std::get<1>(data),
                             std::get<2>(data));
    }
    default:
      throw std::invalid_argument("Invalid Axis::XYZ");
  }
}

}  // namespace xfdtd::transform

#endif  // _XFDTD_CORE_TRANSFORM_H_
