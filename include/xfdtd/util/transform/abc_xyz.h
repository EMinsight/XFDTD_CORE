#ifndef __XFDTD_CORE_TRANSFORM_ABC_XYZ_H__
#define __XFDTD_CORE_TRANSFORM_ABC_XYZ_H__

#include <xfdtd/coordinate_system/coordinate_system.h>

#include <tuple>

namespace xfdtd {

namespace transform {

template <typename T>
inline constexpr auto xYZToABC(const T& x, const T& y, const T& z,
                               Axis::XYZ xyz) -> std::tuple<T, T, T> {
  switch (xyz) {
    case Axis::XYZ::X:
      return {y, z, x};
    case Axis::XYZ::Y:
      return {z, x, y};
    case Axis::XYZ::Z:
      return {x, y, z};
  }
}

template <typename T>
inline constexpr auto xYZToABC(T&& x, T&& y, T&& z,
                               Axis::XYZ xyz) -> std::tuple<T, T, T> {
  switch (xyz) {
    case Axis::XYZ::X:
      return {y, z, x};
    case Axis::XYZ::Y:
      return {z, x, y};
    case Axis::XYZ::Z:
      return {x, y, z};
  }
}

template <typename T>
inline constexpr auto aBCToXYZ(const T& a, const T& b, const T& c,
                               Axis::XYZ xyz) -> std::tuple<T, T, T> {
  switch (xyz) {
    case Axis::XYZ::X:
      return {c, a, b};
    case Axis::XYZ::Y:
      return {b, c, a};
    case Axis::XYZ::Z:
      return {a, b, c};
  }
}

template <typename T>
inline constexpr auto aBCToXYZ(T&& a, T&& b, T&& c,
                               Axis::XYZ xyz) -> std::tuple<T, T, T> {
  switch (xyz) {
    case Axis::XYZ::X:
      return {c, a, b};
    case Axis::XYZ::Y:
      return {b, c, a};
    case Axis::XYZ::Z:
      return {a, b, c};
  }
}

template <typename T, Axis::XYZ xyz>
inline constexpr auto xYZToABC(const T& x, const T& y,
                               const T& z) -> std::tuple<T, T, T> {
  if constexpr (xyz == Axis::XYZ::X) {
    return {y, z, x};
  } else if constexpr (xyz == Axis::XYZ::Y) {
    return {z, x, y};
  } else if constexpr (xyz == Axis::XYZ::Z) {
    return {x, y, z};
  }
}

template <typename T, Axis::XYZ xyz>
inline constexpr auto xYZToABC(T&& x, T&& y, T&& z) -> std::tuple<T, T, T> {
  if constexpr (xyz == Axis::XYZ::X) {
    return {y, z, x};
  } else if constexpr (xyz == Axis::XYZ::Y) {
    return {z, x, y};
  } else if constexpr (xyz == Axis::XYZ::Z) {
    return {x, y, z};
  }
}

template <typename T, Axis::XYZ xyz>
inline constexpr auto aBCToXYZ(const T& a, const T& b,
                               const T& c) -> std::tuple<T, T, T> {
  if constexpr (xyz == Axis::XYZ::X) {
    return {c, a, b};
  } else if constexpr (xyz == Axis::XYZ::Y) {
    return {b, c, a};
  } else if constexpr (xyz == Axis::XYZ::Z) {
    return {a, b, c};
  }
}

template <typename T, Axis::XYZ xyz>
inline constexpr auto aBCToXYZ(T&& a, T&& b, T&& c) -> std::tuple<T, T, T> {
  if constexpr (xyz == Axis::XYZ::X) {
    return {c, a, b};
  } else if constexpr (xyz == Axis::XYZ::Y) {
    return {b, c, a};
  } else if constexpr (xyz == Axis::XYZ::Z) {
    return {a, b, c};
  }
}

template <typename T>
inline auto xYZToABC(const std::tuple<T, T, T>& data, Axis::XYZ xyz) {
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

template <typename T>
inline auto aBCToXYZ(const std::tuple<T, T, T>& data, Axis::XYZ xyz) {
  switch (xyz) {
    case Axis::XYZ::X: {
      // x: c y: a z: b
      return std::make_tuple(std::get<2>(data), std::get<0>(data),
                             std::get<1>(data));
    }
    case Axis::XYZ::Y: {
      // x: b y: c z: a
      return std::make_tuple(std::get<1>(data), std::get<2>(data),
                             std::get<0>(data));
    }
    case Axis::XYZ::Z: {
      // x: a y: b z: c
      return std::make_tuple(std::get<0>(data), std::get<1>(data),
                             std::get<2>(data));
    }
    default:
      throw std::invalid_argument("Invalid Axis::XYZ");
  }
}

}  // namespace transform

}  // namespace xfdtd

#endif  // __XFDTD_CORE_TRANSFORM_ABC_XYZ_H__
