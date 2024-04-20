#ifndef __XFDTD_CORE_TYPE_DEFINE_H__
#define __XFDTD_CORE_TYPE_DEFINE_H__

#include <xtensor/xtensor.hpp>
#include <xtensor/xarray.hpp>
#include <cstddef>

namespace xfdtd {

#if defined XFDTD_CORE_SINGLE_PRECISION
using Real = float;
#else
using Real = double;
#endif

using Index = std::size_t;

template <typename T>
using Array1D = xt::xtensor<T, 1>;

template <typename T>
using Array2D = xt::xtensor<T, 2>;

template <typename T>
using Array3D = xt::xtensor<T, 3>;

template <typename T>
using Array4D = xt::xtensor<T, 4>;

template <typename T>
using Array = xt::xarray<T>;

template <typename T>
concept Number = std::is_arithmetic_v<T>;

template <typename T>
concept ArrLike = requires(T t) {
  { t.size() } -> std::convertible_to<std::size_t>;
  { t[0] } -> Number;
};

namespace unit {

enum class Length {
  Meter,
  Centimeter,
  Millimeter,
  Micrometer,
  Nanometer,
  Angstrom
};

enum class Time { Hour, Minute, Second, Millisecond, Microsecond, Nanosecond };

enum class Frequency { Hertz, Kilohertz, Megahertz, Gigahertz, Terahertz };

}  // namespace unit

}  // namespace xfdtd

#endif  // __XFDTD_CORE_TYPE_DEFINE_H__
