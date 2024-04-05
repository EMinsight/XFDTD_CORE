#include <xfdtd/nffft/nffft.h>
#include <xfdtd/util/constant.h>
#include <xfdtd/util/transform.h>

#include <complex>

#include "nffft/nffft_fd_data.h"
#include "xfdtd/coordinate_system/coordinate_system.h"

namespace xfdtd {

auto FDPlaneData::aTheta(const xt::xtensor<double, 1>& theta,
                         const xt::xtensor<double, 1>& phi,
                         const Vector& origin) const
    -> xt::xtensor<std::complex<double>, 1> {
  return getPotential<Potential::A, transform::SCS::Theta>(theta, phi, origin);
  if (theta.size() == 0 || phi.size() == 0) {
    throw XFDTDNFFFTException("Invalid theta or phi size");
  }

  if (theta.size() != 1 && phi.size() != 1) {
    throw XFDTDNFFFTException("There must be only one theta or phi");
  }

  const auto nun_angle = theta.size() * phi.size();
  auto a_theta = xt::zeros<std::complex<double>>({nun_angle});

  auto&& cos_t{xt::cos(theta)};
  auto&& sin_t{xt::sin(theta)};
  auto&& cos_p{xt::cos(phi)};
  auto&& sin_p{xt::sin(phi)};
  auto&& sin_t_sin_p{sin_t * sin_p};
  auto&& sin_t_cos_p{sin_t * cos_p};
  xt::xtensor<std::complex<double>, 1> a_theta_xn =
      xt::zeros<std::complex<double>>({nun_angle});
  xt::xtensor<std::complex<double>, 1> a_theta_yn =
      xt::zeros<std::complex<double>>({nun_angle});
  xt::zeros<std::complex<double>>({nun_angle});
  xt::xtensor<std::complex<double>, 1> a_theta_zn =
      xt::zeros<std::complex<double>>({nun_angle});
  xt::xtensor<std::complex<double>, 1> a_theta_xp =
      xt::zeros<std::complex<double>>({nun_angle});
  xt::xtensor<std::complex<double>, 1> a_theta_yp =
      xt::zeros<std::complex<double>>({nun_angle});
  xt::xtensor<std::complex<double>, 1> a_theta_zp =
      xt::zeros<std::complex<double>>({nun_angle});

  const auto wave_number = 2.0 * constant::PI * _freq / constant::C_0;

  struct AlwaysOne {
    auto operator()(std::size_t i) const -> auto { return 1; }
  } always_one;

  auto calculate_a_theta =
      [&origin, &sin_t_cos_p, &sin_t_sin_p, &cos_t, &wave_number](
          const auto& task, const auto& node_x, const auto& node_y,
          const auto& node_z, const auto& dx, const auto& dy, const auto& dz,
          const xt::xtensor<double, 1>& transform_a,
          const xt::xtensor<double, 1>& transform_b,
          const xt::xtensor<std::complex<double>, 3>& ja,
          const xt::xtensor<std::complex<double>, 3>& jb, auto& a_theta) {
        if (!task.valid()) {
          return;
        }

        const auto is = task.xRange().start();
        const auto js = task.yRange().start();
        const auto ks = task.zRange().start();
        const auto ie = task.xRange().end();
        const auto je = task.yRange().end();
        const auto ke = task.zRange().end();
        using namespace std::complex_literals;

        for (auto i{is}; i < ie; ++i) {
          for (auto j{js}; j < je; ++j) {
            for (auto k{ks}; k < ke; ++k) {
              auto&& r{Vector{node_x(i), node_y(j), node_z(k)} - origin};
              auto&& ds{dx(i) * dy(j) * dz(k)};
              auto&& phase_shift{xt::exp(
                  1i * wave_number *
                  (r.x() * sin_t_cos_p + r.y() * sin_t_sin_p + r.z() * cos_t))};
              a_theta += (ja(i - is, j - js, k - ks) * transform_a +
                          +jb(i - is, j - js, k - ks) * transform_b) *
                         phase_shift * ds;
            }
          }
        }
      };

  auto&& cos_t_sin_p{cos_t * sin_p};
  auto&& cos_t_cos_p{cos_t * cos_p};
  auto&& n_sin_t{-1 * sin_t};

  calculate_a_theta(_task_xn, _grid_space->eNodeX(), _grid_space->hNodeY(),
                    _grid_space->hNodeZ(), always_one, _grid_space->eSizeY(),
                    _grid_space->eSizeZ(), cos_t_sin_p, n_sin_t, _jy_xn, _jz_xn,
                    a_theta_xn);

  calculate_a_theta(_task_xp, _grid_space->eNodeX(), _grid_space->hNodeY(),
                    _grid_space->hNodeZ(), always_one, _grid_space->eSizeY(),
                    _grid_space->eSizeZ(), cos_t_sin_p, n_sin_t, _jy_xp, _jz_xp,
                    a_theta_xp);

  calculate_a_theta(_task_yn, _grid_space->hNodeX(), _grid_space->eNodeY(),
                    _grid_space->hNodeZ(), _grid_space->eSizeX(), always_one,
                    _grid_space->eSizeZ(), n_sin_t, cos_t_cos_p, _jz_yn, _jx_yn,
                    a_theta_yn);

  calculate_a_theta(_task_yp, _grid_space->hNodeX(), _grid_space->eNodeY(),
                    _grid_space->hNodeZ(), _grid_space->eSizeX(), always_one,
                    _grid_space->eSizeZ(), n_sin_t, cos_t_cos_p, _jz_yp, _jx_yp,
                    a_theta_yp);

  calculate_a_theta(_task_zn, _grid_space->hNodeX(), _grid_space->hNodeY(),
                    _grid_space->eNodeZ(), _grid_space->eSizeX(),
                    _grid_space->eSizeY(), always_one, cos_t_cos_p, cos_t_sin_p,
                    _jx_zn, _jy_zn, a_theta_zn);

  calculate_a_theta(_task_zp, _grid_space->hNodeX(), _grid_space->hNodeY(),
                    _grid_space->eNodeZ(), _grid_space->eSizeX(),
                    _grid_space->eSizeY(), always_one, cos_t_cos_p, cos_t_sin_p,
                    _jx_zp, _jy_zp, a_theta_zp);

  return a_theta_xn + a_theta_xp + a_theta_yn + a_theta_yp + a_theta_zn +
         a_theta_zp;
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
  return calculatePotential<p, Axis::Direction::XN>(
             theta, phi, origin,
             transform::cCSToSCSTransformMatrix<Axis::XYZ::Y, scs>(theta, phi),
             transform::cCSToSCSTransformMatrix<Axis::XYZ::Z, scs>(theta,
                                                                   phi)) +
         calculatePotential<p, Axis::Direction::XP>(
             theta, phi, origin,
             transform::cCSToSCSTransformMatrix<Axis::XYZ::Y, scs>(theta, phi),
             transform::cCSToSCSTransformMatrix<Axis::XYZ::Z, scs>(theta,
                                                                   phi)) +
         calculatePotential<p, Axis::Direction::YN>(
             theta, phi, origin,
             transform::cCSToSCSTransformMatrix<Axis::XYZ::Z, scs>(theta, phi),
             transform::cCSToSCSTransformMatrix<Axis::XYZ::X, scs>(theta,
                                                                   phi)) +
         calculatePotential<p, Axis::Direction::YP>(
             theta, phi, origin,
             transform::cCSToSCSTransformMatrix<Axis::XYZ::Z, scs>(theta, phi),
             transform::cCSToSCSTransformMatrix<Axis::XYZ::X, scs>(theta,
                                                                   phi)) +
         calculatePotential<p, Axis::Direction::ZN>(
             theta, phi, origin,
             transform::cCSToSCSTransformMatrix<Axis::XYZ::X, scs>(theta, phi),
             transform::cCSToSCSTransformMatrix<Axis::XYZ::Y, scs>(theta,
                                                                   phi)) +
         calculatePotential<p, Axis::Direction::ZP>(
             theta, phi, origin,
             transform::cCSToSCSTransformMatrix<Axis::XYZ::X, scs>(theta, phi),
             transform::cCSToSCSTransformMatrix<Axis::XYZ::Y, scs>(theta, phi));
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
    return {};
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
