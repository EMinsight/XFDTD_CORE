#include "nffft/nffft_fd_data.h"

namespace xfdtd {

// Ey: i, j + 1/2, k. Need i, j+ 1/2, k + 1/2. k
static auto interpolateEyFaceX(const auto& ey, const auto& i, const auto& j,
                               const auto& k) {
  return 0.5 * (ey(i, j, k + 1) + ey(i, j, k));
}

// Ez: i, j, k + 1/2. Need i, j + 1/2, k+ 1/2. j
static auto interpolateEzFaceX(const auto& ez, const auto& i, const auto& j,
                               const auto& k) {
  return 0.5 * (ez(i, j + 1, k) + ez(i, j, k));
}

// Hy: i + 1/2, j, k + 1/2. Need i, j + 1/2, k+1/2. i, j
static auto interpolateHyFaceX(const auto& hy, const auto& i, const auto& j,
                               const auto& k) {
  return 0.25 * (hy(i, j + 1, k) + hy(i, j, k) + hy(i - 1, j + 1, k) +
                 hy(i - 1, j, k));
}

// Hz: i + 1/2, j+ 1/2, k. Need i, j + 1/2, k + 1/2. i, k
static auto interpolateHzFaceX(const auto& hz, const auto& i, const auto& j,
                               const auto& k) {
  return 0.25 * (hz(i, j, k + 1) + hz(i, j, k) + hz(i - 1, j, k + 1) +
                 hz(i - 1, j, k));
}

// Ez: i, j, k + 1/2. Need i + 1/2, j, k + 1/2. i
static auto interpolateEzFaceY(const auto& ez, const auto& i, const auto& j,
                               const auto& k) {
  return 0.5 * (ez(i + 1, j, k) + ez(i, j, k));
}

// Ex: i + 1/2, j, k. Need i + 1/2, j, k + 1/2. k
static auto interpolateExFaceY(const auto& ex, const auto& i, const auto& j,
                               const auto& k) {
  return 0.5 * (ex(i, j, k + 1) + ex(i, j, k));
}

// Hz: i + 1/2, j + 1/2, k. Need i + 1/2, j, k + 1/2. j, k
static auto interpolateHzFaceY(const auto& hz, const auto& i, const auto& j,
                               const auto& k) {
  return 0.25 * (hz(i, j, k + 1) + hz(i, j, k) + hz(i, j - 1, k + 1) +
                 +hz(i, j - 1, k));
}

// Hx: i, j + 1/2, k + 1/2. Need i + 1/2, j, k + 1/2. i, j
static auto interpolateHxFaceY(const auto& hx, const auto& i, const auto& j,
                               const auto& k) {
  return 0.25 * (hx(i + 1, j, k) + hx(i, j, k) + hx(i + 1, j - 1, k) +
                 hx(i, j - 1, k));
}

// Ex: i + 1/2, j, k. Need i + 1/2, j + 1/2, k. j
static auto interpolateExFaceZ(const auto& ex, const auto& i, const auto& j,
                               const auto& k) {
  return 0.5 * (ex(i, j + 1, k) + ex(i, j, k));
}

// Ey: i, j + 1/2, k. Need i + 1/2, j + 1/2, k. i
static auto interpolateEyFaceZ(const auto& ey, const auto& i, const auto& j,
                               const auto& k) {
  return 0.5 * (ey(i + 1, j, k) + ey(i, j, k));
}

// Hx: i, j + 1/2, k + 1/2. Need i + 1/2, j + 1/2, k. i, k
static auto interpolateHxFaceZ(const auto& hx, const auto& i, const auto& j,
                               const auto& k) {
  return 0.25 * (hx(i + 1, j, k) + hx(i, j, k) + hx(i, j, k - 1) +
                 hx(i + 1, j, k - 1));
}

// Hy: i + 1/2, j, k + 1/2. Need i + 1/2, j + 1/2, k. j, k
static auto interpolateHyFaceZ(const auto& hy, const auto& i, const auto& j,
                               const auto& k) {
  return 0.25 * (hy(i, j, k) + hy(i, j + 1, k) + hy(i, j, k - 1) +
                 hy(i, j + 1, k - 1));
}

auto FDPlaneData::update(std::size_t current_time_step) -> void {
  calculateJ<Axis::Direction::XN>(current_time_step);
  calculateJ<Axis::Direction::XP>(current_time_step);
  calculateJ<Axis::Direction::YN>(current_time_step);
  calculateJ<Axis::Direction::YP>(current_time_step);
  calculateJ<Axis::Direction::ZN>(current_time_step);
  calculateJ<Axis::Direction::ZP>(current_time_step);

  calculateM<Axis::Direction::XN>(current_time_step);
  calculateM<Axis::Direction::XP>(current_time_step);
  calculateM<Axis::Direction::YN>(current_time_step);
  calculateM<Axis::Direction::YP>(current_time_step);
  calculateM<Axis::Direction::ZN>(current_time_step);
  calculateM<Axis::Direction::ZP>(current_time_step);
}

template <Axis::XYZ xyz, EMF::Field f>
auto FDPlaneData::interpolate(const std::size_t& i, const std::size_t& j,
                              const std::size_t& k) -> Real {
  if constexpr (xyz == Axis::XYZ::X) {
    if constexpr (f == EMF::Field::HZ) {
      return interpolateHzFaceX(_emf->hz(), i, j, k);
    } else if constexpr (f == EMF::Field::HY) {
      return interpolateHyFaceX(_emf->hy(), i, j, k);
    } else if constexpr (f == EMF::Field::EZ) {
      return interpolateEzFaceX(_emf->ez(), i, j, k);
    } else if constexpr (f == EMF::Field::EY) {
      return interpolateEyFaceX(_emf->ey(), i, j, k);
    } else {
      throw XFDTDNFFFTException("Invalid field");
    }
  } else if constexpr (xyz == Axis::XYZ::Y) {
    if constexpr (f == EMF::Field::HZ) {
      return interpolateHzFaceY(_emf->hz(), i, j, k);
    } else if constexpr (f == EMF::Field::HX) {
      return interpolateHxFaceY(_emf->hx(), i, j, k);
    } else if constexpr (f == EMF::Field::EZ) {
      return interpolateEzFaceY(_emf->ez(), i, j, k);
    } else if constexpr (f == EMF::Field::EX) {
      return interpolateExFaceY(_emf->ex(), i, j, k);
    } else {
      throw XFDTDNFFFTException("Invalid field");
    }
  } else if constexpr (xyz == Axis::XYZ::Z) {
    if constexpr (f == EMF::Field::HX) {
      return interpolateHxFaceZ(_emf->hx(), i, j, k);
    } else if constexpr (f == EMF::Field::HY) {
      return interpolateHyFaceZ(_emf->hy(), i, j, k);
    } else if constexpr (f == EMF::Field::EX) {
      return interpolateExFaceZ(_emf->ex(), i, j, k);
    } else if constexpr (f == EMF::Field::EY) {
      return interpolateEyFaceZ(_emf->ey(), i, j, k);
    } else {
      throw XFDTDNFFFTException("Invalid field");
    }
  } else {
    throw XFDTDNFFFTException("Invalid field");
  }
}

template <Axis::Direction direction>
auto FDPlaneData::calculateM(std::size_t current_time_step) -> void {
  const auto task = this->task<direction>();

  if (!task.valid()) {
    return;
  }

  Real coff_ma = (Axis::directionPositive<direction>()) ? 1.0 : -1.0;
  Real coff_mb = (Axis::directionPositive<direction>()) ? -1.0 : 1.0;

  const auto is = task.xRange().start();
  const auto js = task.yRange().start();
  const auto ks = task.zRange().start();
  const auto ie = task.xRange().end();
  const auto je = task.yRange().end();
  const auto ke = task.zRange().end();

  auto [ma, mb] = surfaceM<direction>();

  constexpr auto xyz = Axis::fromDirectionToXYZ<direction>();

  constexpr auto ea = tangentialEaFieldEnum<xyz>();
  constexpr auto eb = tangentialEbFieldEnum<xyz>();

  for (auto i{is}; i < ie; ++i) {
    for (auto j{js}; j < je; ++j) {
      for (auto k{ks}; k < ke; ++k) {
        ma(i - is, j - js, k - ks) += coff_ma * interpolate<xyz, eb>(i, j, k) *
                                      _transform_e(current_time_step);
        mb(i - is, j - js, k - ks) += coff_mb * interpolate<xyz, ea>(i, j, k) *
                                      _transform_e(current_time_step);
      }
    }
  }
}

template <Axis::Direction direction>
auto FDPlaneData::calculateJ(std::size_t current_time_step) -> void {
  const auto task = this->task<direction>();

  if (!task.valid()) {
    return;
  }

  constexpr Real coff_ja =
      (Axis::directionPositive<direction>()) ? -1.0 : 1.0;
  constexpr Real coff_jb =
      (Axis::directionPositive<direction>()) ? 1.0 : -1.0;

  constexpr auto xyz = Axis::fromDirectionToXYZ<direction>();

  const auto is = task.xRange().start();
  const auto js = task.yRange().start();
  const auto ks = task.zRange().start();
  const auto ie = task.xRange().end();
  const auto je = task.yRange().end();
  const auto ke = task.zRange().end();

  auto [ja, jb] = surfaceJ<direction>();

  constexpr auto ha = tangentialHaFieldEnum<xyz>();
  constexpr auto hb = tangentialHbFieldEnum<xyz>();

  for (auto i{is}; i < ie; ++i) {
    for (auto j{js}; j < je; ++j) {
      for (auto k{ks}; k < ke; ++k) {
        ja(i - is, j - js, k - ks) += coff_ja * interpolate<xyz, hb>(i, j, k) *
                                      _transform_h(current_time_step);
        jb(i - is, j - js, k - ks) += coff_jb * interpolate<xyz, ha>(i, j, k) *
                                      _transform_h(current_time_step);
      }
    }
  }
}

}  // namespace xfdtd
