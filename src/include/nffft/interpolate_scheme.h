#ifndef __XFDTD_CORE_INTERPOLATE_SCHEME_H__
#define __XFDTD_CORE_INTERPOLATE_SCHEME_H__

#include <utility>

#include "xfdtd/common/type_define.h"
#include "xfdtd/electromagnetic_field/electromagnetic_field.h"
#include "xfdtd/exception/exception.h"

namespace xfdtd::interpolate {

class XFDTDInterpolateException : public XFDTDException {
 public:
  explicit XFDTDInterpolateException(std::string message)
      : XFDTDException{std::move(message)} {};
};

// Ey: i, j + 1/2, k. Need i, j+ 1/2, k + 1/2. k
inline auto interpolateEyFaceX(const auto& ey, const auto& i, const auto& j,
                               const auto& k) {
  return 0.5 * (ey(i, j, k + 1) + ey(i, j, k));
}

// Ez: i, j, k + 1/2. Need i, j + 1/2, k+ 1/2. j
inline auto interpolateEzFaceX(const auto& ez, const auto& i, const auto& j,
                               const auto& k) {
  return 0.5 * (ez(i, j + 1, k) + ez(i, j, k));
}

// Hy: i + 1/2, j, k + 1/2. Need i, j + 1/2, k+1/2. i, j
inline auto interpolateHyFaceX(const auto& hy, const auto& i, const auto& j,
                               const auto& k) {
  return 0.25 * (hy(i, j + 1, k) + hy(i, j, k) + hy(i - 1, j + 1, k) +
                 hy(i - 1, j, k));
}

// Hz: i + 1/2, j+ 1/2, k. Need i, j + 1/2, k + 1/2. i, k
inline auto interpolateHzFaceX(const auto& hz, const auto& i, const auto& j,
                               const auto& k) {
  return 0.25 * (hz(i, j, k + 1) + hz(i, j, k) + hz(i - 1, j, k + 1) +
                 hz(i - 1, j, k));
}

// Ez: i, j, k + 1/2. Need i + 1/2, j, k + 1/2. i
inline auto interpolateEzFaceY(const auto& ez, const auto& i, const auto& j,
                               const auto& k) {
  return 0.5 * (ez(i + 1, j, k) + ez(i, j, k));
}

// Ex: i + 1/2, j, k. Need i + 1/2, j, k + 1/2. k
inline auto interpolateExFaceY(const auto& ex, const auto& i, const auto& j,
                               const auto& k) {
  return 0.5 * (ex(i, j, k + 1) + ex(i, j, k));
}

// Hz: i + 1/2, j + 1/2, k. Need i + 1/2, j, k + 1/2. j, k
inline auto interpolateHzFaceY(const auto& hz, const auto& i, const auto& j,
                               const auto& k) {
  return 0.25 * (hz(i, j, k + 1) + hz(i, j, k) + hz(i, j - 1, k + 1) +
                 +hz(i, j - 1, k));
}

// Hx: i, j + 1/2, k + 1/2. Need i + 1/2, j, k + 1/2. i, j
inline auto interpolateHxFaceY(const auto& hx, const auto& i, const auto& j,
                               const auto& k) {
  return 0.25 * (hx(i + 1, j, k) + hx(i, j, k) + hx(i + 1, j - 1, k) +
                 hx(i, j - 1, k));
}

// Ex: i + 1/2, j, k. Need i + 1/2, j + 1/2, k. j
inline auto interpolateExFaceZ(const auto& ex, const auto& i, const auto& j,
                               const auto& k) {
  return 0.5 * (ex(i, j + 1, k) + ex(i, j, k));
}

// Ey: i, j + 1/2, k. Need i + 1/2, j + 1/2, k. i
inline auto interpolateEyFaceZ(const auto& ey, const auto& i, const auto& j,
                               const auto& k) {
  return 0.5 * (ey(i + 1, j, k) + ey(i, j, k));
}

// Hx: i, j + 1/2, k + 1/2. Need i + 1/2, j + 1/2, k. i, k
inline auto interpolateHxFaceZ(const auto& hx, const auto& i, const auto& j,
                               const auto& k) {
  return 0.25 * (hx(i + 1, j, k) + hx(i, j, k) + hx(i, j, k - 1) +
                 hx(i + 1, j, k - 1));
}

// Hy: i + 1/2, j, k + 1/2. Need i + 1/2, j + 1/2, k. j, k
inline auto interpolateHyFaceZ(const auto& hy, const auto& i, const auto& j,
                               const auto& k) {
  return 0.25 * (hy(i, j, k) + hy(i, j + 1, k) + hy(i, j, k - 1) +
                 hy(i, j + 1, k - 1));
}

template <Axis::XYZ xyz, EMF::Field f>
inline auto interpolateSurfaceCenter(const auto& data, Index i, Index j,
                                     Index k) {
  if constexpr (xyz == Axis::XYZ::X) {
    if constexpr (f == EMF::Field::HZ) {
      return interpolateHzFaceX(data, i, j, k);
    } else if constexpr (f == EMF::Field::HY) {
      return interpolateHyFaceX(data, i, j, k);
    } else if constexpr (f == EMF::Field::EZ) {
      return interpolateEzFaceX(data, i, j, k);
    } else if constexpr (f == EMF::Field::EY) {
      return interpolateEyFaceX(data, i, j, k);
    } else {
      throw XFDTDInterpolateException("Invalid field");
    }
  } else if constexpr (xyz == Axis::XYZ::Y) {
    if constexpr (f == EMF::Field::HZ) {
      return interpolateHzFaceY(data, i, j, k);
    } else if constexpr (f == EMF::Field::HX) {
      return interpolateHxFaceY(data, i, j, k);
    } else if constexpr (f == EMF::Field::EZ) {
      return interpolateEzFaceY(data, i, j, k);
    } else if constexpr (f == EMF::Field::EX) {
      return interpolateExFaceY(data, i, j, k);
    } else {
      throw XFDTDInterpolateException("Invalid field");
    }
  } else if constexpr (xyz == Axis::XYZ::Z) {
    if constexpr (f == EMF::Field::HX) {
      return interpolateHxFaceZ(data, i, j, k);
    } else if constexpr (f == EMF::Field::HY) {
      return interpolateHyFaceZ(data, i, j, k);
    } else if constexpr (f == EMF::Field::EX) {
      return interpolateExFaceZ(data, i, j, k);
    } else if constexpr (f == EMF::Field::EY) {
      return interpolateEyFaceZ(data, i, j, k);
    } else {
      throw XFDTDInterpolateException("Invalid field");
    }
  } else {
    throw XFDTDInterpolateException("Invalid field");
  }
}

// template <Axis::XYZ xyz, EMF::Field f>
// inline auto interpolate(const EMF* emf, Index i, Index j, Index k) {
//   if constexpr (xyz == Axis::XYZ::X) {
//     if constexpr (f == EMF::Field::HZ) {
//       return interpolateHzFaceX(emf->hz(), i, j, k);
//     } else if constexpr (f == EMF::Field::HY) {
//       return interpolateHyFaceX(emf->hy(), i, j, k);
//     } else if constexpr (f == EMF::Field::EZ) {
//       return interpolateEzFaceX(emf->ez(), i, j, k);
//     } else if constexpr (f == EMF::Field::EY) {
//       return interpolateEyFaceX(emf->ey(), i, j, k);
//     } else {
//       throw XFDTDInterpolateException("Invalid field");
//     }
//   } else if constexpr (xyz == Axis::XYZ::Y) {
//     if constexpr (f == EMF::Field::HZ) {
//       return interpolateHzFaceY(emf->hz(), i, j, k);
//     } else if constexpr (f == EMF::Field::HX) {
//       return interpolateHxFaceY(emf->hx(), i, j, k);
//     } else if constexpr (f == EMF::Field::EZ) {
//       return interpolateEzFaceY(emf->ez(), i, j, k);
//     } else if constexpr (f == EMF::Field::EX) {
//       return interpolateExFaceY(emf->ex(), i, j, k);
//     } else {
//       throw XFDTDInterpolateException("Invalid field");
//     }
//   } else if constexpr (xyz == Axis::XYZ::Z) {
//     if constexpr (f == EMF::Field::HX) {
//       return interpolateHxFaceZ(emf->hx(), i, j, k);
//     } else if constexpr (f == EMF::Field::HY) {
//       return interpolateHyFaceZ(emf->hy(), i, j, k);
//     } else if constexpr (f == EMF::Field::EX) {
//       return interpolateExFaceZ(emf->ex(), i, j, k);
//     } else if constexpr (f == EMF::Field::EY) {
//       return interpolateEyFaceZ(emf->ey(), i, j, k);
//     } else {
//       throw XFDTDInterpolateException("Invalid field");
//     }
//   } else {
//     throw XFDTDInterpolateException("Invalid field");
//   }
// }

}  // namespace xfdtd::interpolate

#endif  // __XFDTD_CORE_INTERPOLATE_SCHEME_H__
