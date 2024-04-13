#ifndef _XFDTD_CORE_CONSTANT_H_
#define _XFDTD_CORE_CONSTANT_H_

#include <xfdtd/common/type_define.h>

#include <complex>
#include <numbers>

namespace xfdtd::constant {

inline constexpr Real PI{std::numbers::pi_v<Real>};
inline constexpr Real EPSILON_0{8.854187817e-12};
inline constexpr Real MU_0{4 * PI * 1e-7};
inline const Real C_0{1 / std::sqrt(EPSILON_0 * MU_0)};
inline const Real Z_0{std::sqrt(MU_0 / EPSILON_0)};

inline constexpr Real SIGMA_E_ZERO_APPROX{1e-20};
inline constexpr Real SIGMA_M_ZERO_APPROX{1e-20};

inline constexpr Real INF{std::numeric_limits<Real>::infinity()};
inline constexpr Real NEG_INF{-1 * std::numeric_limits<Real>::infinity()};

constexpr auto II = std::complex<Real>(0.0, 1.0);

}  // namespace xfdtd::constant

#endif  // _XFDTD_CORE_CONSTANT_H_
