#ifndef _XFDTD_LIB_CONSTANT_H_
#define _XFDTD_LIB_CONSTANT_H_

#ifdef _MSC_VER
#define _USE_MATH_DEFINES
#endif

#include <cmath>

namespace xfdtd::constant {

inline constexpr double PI{M_PI};
inline constexpr double EPSILON_0{8.854187817e-12};
inline constexpr double MU_0{4 * PI * 1e-7};
inline const double C_0{1 / std::sqrt(EPSILON_0 * MU_0)};
inline const double Z_0{std::sqrt(MU_0 / EPSILON_0)};

}  // namespace xfdtd::constant

#endif  // _XFDTD_LIB_CONSTANT_H_
