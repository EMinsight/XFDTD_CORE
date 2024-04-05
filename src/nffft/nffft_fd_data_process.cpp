#include <xfdtd/coordinate_system/coordinate_system.h>
#include <xfdtd/nffft/nffft.h>
#include <xfdtd/util/constant.h>
#include <xfdtd/util/transform.h>

#include <complex>
#include <future>
#include <vector>

#include "nffft/nffft_fd_data.h"

namespace xfdtd {

auto FDPlaneData::aTheta(const xt::xtensor<double, 1>& theta,
                         const xt::xtensor<double, 1>& phi,
                         const Vector& origin) const
    -> xt::xtensor<std::complex<double>, 1> {
  return getPotential<Potential::A, transform::SCS::Theta>(theta, phi, origin);
}

auto FDPlaneData::fPhi(const xt::xtensor<double, 1>& theta,
                       const xt::xtensor<double, 1>& phi,
                       const Vector& origin) const
    -> xt::xtensor<std::complex<double>, 1> {
  return getPotential<Potential::F, transform::SCS::Phi>(theta, phi, origin);
}

auto FDPlaneData::aPhi(const xt::xtensor<double, 1>& theta,
                       const xt::xtensor<double, 1>& phi,
                       const Vector& origin) const
    -> xt::xtensor<std::complex<double>, 1> {
  return getPotential<Potential::A, transform::SCS::Phi>(theta, phi, origin);
}

auto FDPlaneData::fTheta(const xt::xtensor<double, 1>& theta,
                         const xt::xtensor<double, 1>& phi,
                         const Vector& origin) const
    -> xt::xtensor<std::complex<double>, 1> {
  return getPotential<Potential::F, transform::SCS::Theta>(theta, phi, origin);
}

auto FDPlaneData::power() const -> double {
  std::vector<std::future<std::complex<double>>> res;
  res.emplace_back(
      std::async(&FDPlaneData::calculatePower<Axis::Direction::XN>, this));
  res.emplace_back(
      std::async(&FDPlaneData::calculatePower<Axis::Direction::XP>, this));
  res.emplace_back(
      std::async(&FDPlaneData::calculatePower<Axis::Direction::YN>, this));
  res.emplace_back(
      std::async(&FDPlaneData::calculatePower<Axis::Direction::YP>, this));
  res.emplace_back(
      std::async(&FDPlaneData::calculatePower<Axis::Direction::ZN>, this));
  res.emplace_back(
      std::async(&FDPlaneData::calculatePower<Axis::Direction::ZP>, this));

  return 0.5 * std::real(std::accumulate(
                   res.begin(), res.end(), std::complex<double>{0.0},
                   [](const auto& a, auto&& b) { return a + b.get(); }));

  return 0.5 * std::real(calculatePower<Axis::Direction::XN>() +
                         calculatePower<Axis::Direction::XP>() +
                         calculatePower<Axis::Direction::YN>() +
                         calculatePower<Axis::Direction::YP>() +
                         calculatePower<Axis::Direction::ZN>() +
                         calculatePower<Axis::Direction::ZP>());
}

template <FDPlaneData::Potential p, transform::SCS scs>
auto FDPlaneData::getPotential(const xt::xtensor<double, 1>& theta,
                               const xt::xtensor<double, 1>& phi,
                               const Vector& origin) const
    -> xt::xtensor<std::complex<double>, 1> {
  auto xn = std::async(
      &FDPlaneData::calculatePotential<p, Axis::Direction::XN>, this, theta,
      phi, origin,
      transform::cCSToSCSTransformMatrix<Axis::XYZ::Y, scs>(theta, phi),
      transform::cCSToSCSTransformMatrix<Axis::XYZ::Z, scs>(theta, phi));
  auto xp = std::async(
      &FDPlaneData::calculatePotential<p, Axis::Direction::XP>, this, theta,
      phi, origin,
      transform::cCSToSCSTransformMatrix<Axis::XYZ::Y, scs>(theta, phi),
      transform::cCSToSCSTransformMatrix<Axis::XYZ::Z, scs>(theta, phi));
  auto yn = std::async(
      &FDPlaneData::calculatePotential<p, Axis::Direction::YN>, this, theta,
      phi, origin,
      transform::cCSToSCSTransformMatrix<Axis::XYZ::Z, scs>(theta, phi),
      transform::cCSToSCSTransformMatrix<Axis::XYZ::X, scs>(theta, phi));
  auto yp = std::async(
      &FDPlaneData::calculatePotential<p, Axis::Direction::YP>, this, theta,
      phi, origin,
      transform::cCSToSCSTransformMatrix<Axis::XYZ::Z, scs>(theta, phi),
      transform::cCSToSCSTransformMatrix<Axis::XYZ::X, scs>(theta, phi));
  auto zn = std::async(
      &FDPlaneData::calculatePotential<p, Axis::Direction::ZN>, this, theta,
      phi, origin,
      transform::cCSToSCSTransformMatrix<Axis::XYZ::X, scs>(theta, phi),
      transform::cCSToSCSTransformMatrix<Axis::XYZ::Y, scs>(theta, phi));
  auto zp = std::async(
      &FDPlaneData::calculatePotential<p, Axis::Direction::ZP>, this, theta,
      phi, origin,
      transform::cCSToSCSTransformMatrix<Axis::XYZ::X, scs>(theta, phi),
      transform::cCSToSCSTransformMatrix<Axis::XYZ::Y, scs>(theta, phi));

  return xn.get() + xp.get() + yn.get() + yp.get() + zn.get() + zp.get();
}

template <FDPlaneData::Potential potential, Axis::Direction direction>
auto FDPlaneData::calculatePotential(
    const xt::xtensor<double, 1>& theta, const xt::xtensor<double, 1>& phi,
    const Vector& origin,
    const xt::xtensor<std::complex<double>, 1>& transform_a,
    const xt::xtensor<std::complex<double>, 1>& transform_b) const
    -> xt::xtensor<std::complex<double>, 1> {
  if (theta.size() == 0 || phi.size() == 0) {
    throw XFDTDNFFFTException("Invalid theta or phi size");
  }

  if (theta.size() != 1 && phi.size() != 1) {
    throw XFDTDNFFFTException("There must be only one theta or phi");
  }

  const auto num_angle = theta.size() * phi.size();
  xt::xtensor<std::complex<double>, 1> data =
      xt::zeros<std::complex<double>>({num_angle});

  const auto& task = this->task<direction>();
  if (!task.valid()) {
    return data;
  }

  const auto is = task.xRange().start();
  const auto js = task.yRange().start();
  const auto ks = task.zRange().start();
  const auto ie = task.xRange().end();
  const auto je = task.yRange().end();
  const auto ke = task.zRange().end();

  using namespace std::complex_literals;

  auto wave_number = 2.0 * constant::PI * _freq / constant::C_0;
  auto&& cos_t = xt::cos(theta);
  auto&& sin_t{xt::sin(theta)};
  auto&& cos_p{xt::cos(phi)};
  auto&& sin_p{xt::sin(phi)};
  auto&& sin_t_sin_p{sin_t * sin_p};
  auto&& sin_t_cos_p{sin_t * cos_p};

  constexpr auto xyz = Axis::fromDirectionToXYZ<direction>();

  auto [a, b] = surfaceCurrent<potential, direction>();

  for (auto i{is}; i < ie; ++i) {
    for (auto j{js}; j < je; ++j) {
      for (auto k{ks}; k < ke; ++k) {
        auto&& r = rVector<xyz>(i, j, k) - origin;
        auto&& ds = this->ds<xyz>(i, j, k);
        auto&& phase_shift = xt::exp(
            1i * wave_number *
            (r.x() * sin_t_cos_p + r.y() * sin_t_sin_p + r.z() * cos_t));
        data += (a(i - is, j - js, k - ks) * transform_a +
                 +b(i - is, j - js, k - ks) * transform_b) *
                phase_shift * ds;
      }
    }
  }

  return data;
}

template <Axis::Direction direction>
auto FDPlaneData::calculatePower() const -> std::complex<double> {
  const auto& task = this->task<direction>();
  if (!task.valid()) {
    return 0.0;
  }

  const auto is = task.xRange().start();
  const auto js = task.yRange().start();
  const auto ks = task.zRange().start();
  const auto ie = task.xRange().end();
  const auto je = task.yRange().end();
  const auto ke = task.zRange().end();

  auto power = std::complex<double>{0.0};

  auto [ja, jb] = surfaceJ<direction>();
  auto [ma, mb] = surfaceM<direction>();

  // Power = /int J^* \times M  \cdot ds
  for (auto i{is}; i < ie; ++i) {
    for (auto j{js}; j < je; ++j) {
      for (auto k{ks}; k < ke; ++k) {
        const auto& ds =
            this->ds<Axis::fromDirectionToXYZ<direction>()>(i, j, k);
        power += ds * (std::conj(ja(i - is, j - js, k - ks)) *
                           mb(i - is, j - js, k - ks) -
                       std::conj(jb(i - is, j - js, k - ks)) *
                           ma(i - is, j - js, k - ks));
      }
    }
  }

  constexpr double negative = Axis::directionNegative<direction>() ? -1.0 : 1.0;
  return negative * power;
}

template <Axis::XYZ xyz>
auto FDPlaneData::rVector(size_t i, std::size_t j, std::size_t k) const
    -> Vector {
  if constexpr (xyz == Axis::XYZ::X) {
    return Vector{_grid_space->eNodeX()(i), _grid_space->hNodeY()(j),
                  _grid_space->hNodeZ()(k)};

  } else if constexpr (xyz == Axis::XYZ::Y) {
    return Vector{_grid_space->hNodeX()(i), _grid_space->eNodeY()(j),
                  _grid_space->hNodeZ()(k)};
  } else if constexpr (xyz == Axis::XYZ::Z) {
    return Vector{_grid_space->hNodeX()(i), _grid_space->hNodeY()(j),
                  _grid_space->eNodeZ()(k)};
  } else {
    throw XFDTDNFFFTException("Invalid xyz");
  }
}

template <Axis::XYZ xyz>
auto FDPlaneData::ds(std::size_t i, std::size_t j, std::size_t k) const
    -> double {
  if constexpr (xyz == Axis::XYZ::X) {
    return _grid_space->eSizeY()(j) * _grid_space->eSizeZ()(k);
  } else if constexpr (xyz == Axis::XYZ::Y) {
    return _grid_space->eSizeZ()(k) * _grid_space->eSizeX()(i);
  } else if constexpr (xyz == Axis::XYZ::Z) {
    return _grid_space->eSizeX()(i) * _grid_space->eSizeY()(j);
  } else {
    throw XFDTDNFFFTException("Invalid xyz");
  }
}

}  // namespace xfdtd
