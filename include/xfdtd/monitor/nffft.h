#ifndef _XFDTD_LIB_NFFFT_H_
#define _XFDTD_LIB_NFFFT_H_

#include <complex>
#include <cstddef>
#include <cstdlib>
#include <filesystem>
#include <utility>
#include <xtensor/xnpy.hpp>

#include "xfdtd/calculation_param/calculation_param.h"
#include "xfdtd/coordinate_system/coordinate_system.h"
#include "xfdtd/electromagnetic_field/electromagnetic_field.h"
#include "xfdtd/grid_space/grid_space.h"
#include "xfdtd/util/constant.h"

namespace xfdtd {

class NFFFT {
 public:
  NFFFT(std::size_t distance_x, std::size_t distance_y, std::size_t distance_z,
        xt::xarray<double> frequencies, std::string output_dir);

  NFFFT(const NFFFT&) = delete;

  NFFFT(NFFFT&&) noexcept = default;

  NFFFT& operator=(const NFFFT&) = delete;

  NFFFT& operator=(NFFFT&&) noexcept = default;

  ~NFFFT() = default;

  void init(std::shared_ptr<const GridSpace> grid_space,
            std::shared_ptr<const CalculationParam> calculation_param,
            std::shared_ptr<const EMF> emf);

  void update();

  void output();

  void processFarField(const xt::xarray<double>& theta,
                       const xt::xarray<double>& phi, const Vector& origin);

 protected:
  const GridSpace* gridSpacePtr() const;

  const CalculationParam* calculationParamPtr() const;

  const EMF* emfPtr() const;

  const GridBox* gridBoxPtr() const;

 private:
  std::size_t _distance_x, _distance_y, _distance_z;
  xt::xarray<double> _frequencies;
  std::string _output_dir;

  std::shared_ptr<const GridSpace> _grid_space;
  std::shared_ptr<const CalculationParam> _calculation_param;
  std::shared_ptr<const EMF> _emf;

  std::unique_ptr<GridBox> _grid_box;

  xt::xarray<std::complex<double>> _jx_yn, _jx_yp, _jx_zn, _jx_zp;
  xt::xarray<std::complex<double>> _jy_xn, _jy_xp, _jy_zn, _jy_zp;
  xt::xarray<std::complex<double>> _jz_xn, _jz_xp, _jz_yn, _jz_yp;
  xt::xarray<std::complex<double>> _mx_yn, _mx_yp, _mx_zn, _mx_zp;
  xt::xarray<std::complex<double>> _my_xn, _my_xp, _my_zn, _my_zp;
  xt::xarray<std::complex<double>> _mz_xn, _mz_xp, _mz_yn, _mz_yp;

  xt::xarray<std::complex<double>> _transform_e, _transform_h;

  std::complex<double> _a_x, _a_y, _a_z, _f_x, _f_y, _f_z;
  xt::xarray<std::complex<double>> _a_theta, _a_phi, _f_theta, _f_phi;

  void initDFT();
};

void NFFFT::init(std::shared_ptr<const GridSpace> grid_space,
                 std::shared_ptr<const CalculationParam> calculation_param,
                 std::shared_ptr<const EMF> emf) {
  _grid_space = std::move(grid_space);
  _calculation_param = std::move(calculation_param);
  _emf = std::move(emf);

  auto nx{_grid_space->sizeX()};
  auto ny{_grid_space->sizeY()};
  auto nz{_grid_space->sizeZ()};

  auto is{_distance_x};
  auto js{_distance_y};
  auto ks{_distance_z};

  auto ie{nx - _distance_x};
  auto je{ny - _distance_y};
  auto ke{nz - _distance_z};

  auto size_i{ie - is};
  auto size_j{je - js};
  auto size_k{ke - ks};

  _grid_box =
      std::make_unique<GridBox>(Grid(is, js, ks), Grid(size_i, size_j, size_k));

  _jx_yn = xt::zeros<std::complex<double>>(
      {_frequencies.size(), size_i, 1UL, size_k});
  _jx_yp = xt::zeros<std::complex<double>>(
      {_frequencies.size(), size_i, 1UL, size_k});
  _jx_zn = xt::zeros<std::complex<double>>(
      {_frequencies.size(), size_i, size_j, 1UL});
  _jx_zp = xt::zeros<std::complex<double>>(
      {_frequencies.size(), size_i, size_j, 1UL});
  _jy_xn = xt::zeros<std::complex<double>>(
      {_frequencies.size(), 1UL, size_j, size_k});
  _jy_xp = xt::zeros<std::complex<double>>(
      {_frequencies.size(), 1UL, size_j, size_k});
  _jy_zn = xt::zeros<std::complex<double>>(
      {_frequencies.size(), size_i, size_j, 1UL});
  _jy_zp = xt::zeros<std::complex<double>>(
      {_frequencies.size(), size_i, size_j, 1UL});
  _jz_xn = xt::zeros<std::complex<double>>(
      {_frequencies.size(), 1UL, size_j, size_k});
  _jz_xp = xt::zeros<std::complex<double>>(
      {_frequencies.size(), 1UL, size_j, size_k});
  _jz_yn = xt::zeros<std::complex<double>>(
      {_frequencies.size(), size_i, 1UL, size_k});
  _jz_yp = xt::zeros<std::complex<double>>(
      {_frequencies.size(), size_i, 1UL, size_k});
  _mx_yn = xt::zeros<std::complex<double>>(
      {_frequencies.size(), size_i, 1UL, size_k});
  _mx_yp = xt::zeros<std::complex<double>>(
      {_frequencies.size(), size_i, 1UL, size_k});
  _mx_zn = xt::zeros<std::complex<double>>(
      {_frequencies.size(), size_i, size_j, 1UL});
  _mx_zp = xt::zeros<std::complex<double>>(
      {_frequencies.size(), size_i, size_j, 1UL});
  _my_xn = xt::zeros<std::complex<double>>(
      {_frequencies.size(), 1UL, size_j, size_k});
  _my_xp = xt::zeros<std::complex<double>>(
      {_frequencies.size(), 1UL, size_j, size_k});
  _my_zn = xt::zeros<std::complex<double>>(
      {_frequencies.size(), size_i, size_j, 1UL});
  _my_zp = xt::zeros<std::complex<double>>(
      {_frequencies.size(), size_i, size_j, 1UL});
  _mz_xn = xt::zeros<std::complex<double>>(
      {_frequencies.size(), 1UL, size_j, size_k});
  _mz_xp = xt::zeros<std::complex<double>>(
      {_frequencies.size(), 1UL, size_j, size_k});
  _mz_yn = xt::zeros<std::complex<double>>(
      {_frequencies.size(), size_i, 1UL, size_k});
  _mz_yp = xt::zeros<std::complex<double>>(
      {_frequencies.size(), size_i, 1UL, size_k});
  initDFT();
}

void NFFFT::initDFT() {
  using namespace std::complex_literals;
  const auto total_time_step{
      calculationParamPtr()->timeParam()->endTimeStep() -
      calculationParamPtr()->timeParam()->startTimeStep() + 1};
  const auto dt{calculationParamPtr()->timeParam()->dt()};
  _transform_e =
      xt::zeros<std::complex<double>>({_frequencies.size(), total_time_step});
  _transform_h =
      xt::zeros<std::complex<double>>({_frequencies.size(), total_time_step});
  for (auto f{0}; f < _frequencies.size(); ++f) {
    for (std::size_t t{0}; t < total_time_step; ++t) {
      _transform_e(f, t) =
          dt * std::exp(-1i * 2.0 * constant::PI * (_frequencies(f) * t * dt));
      _transform_h(f, t) = dt * std::exp(-1i * 2.0 * constant::PI *
                                         (_frequencies(f) * (t - 0.5) * dt));
    }
  }
}

void NFFFT::update() {
  auto emf{emfPtr()};
  auto box{gridBoxPtr()};
  auto is{box->origin().i()};
  auto js{box->origin().j()};
  auto ks{box->origin().k()};
  auto ie{box->end().i()};
  auto je{box->end().j()};
  auto ke{box->end().k()};
  auto current_time_step{calculationParamPtr()->timeParam()->currentTimeStep()};

  for (auto j{js}; j < je; ++j) {
    for (auto k{ks}; k < ke; ++k) {
      auto my_xn{-0.5 * (emf->ez()(is, j + 1, k) + emf->ez()(is, j, k))};
      auto mz_xn{0.5 * (emf->ey()(is, j, k + 1) + emf->ey()(is, j, k))};
      auto jy_xn{0.25 *
                 (emf->hz()(is, j, k + 1) + emf->hz()(is, j, k) +
                  emf->hz()(is - 1, j, k + 1) + emf->hz()(is - 1, j, k))};
      auto jz_xn{-0.25 *
                 (emf->hy()(is, j + 1, k) + emf->hy()(is, j, k) +
                  emf->hy()(is - 1, j + 1, k) + emf->hy()(is - 1, j, k))};

      auto my_xp{emf->ezFaceCenterX(ie, j, k)};
      auto mz_xp{-1.0 * emf->eyFaceCenterX(ie, j, k)};
      auto jy_xp{-1.0 * emf->hzFaceCenterX(ie, j, k)};
      auto jz_xp{emf->hyFaceCenterX(ie, j, k)};

      for (int f{0}; f < _frequencies.size(); ++f) {
        _my_xn(f, 0, j - js, k - ks) +=
            my_xn * _transform_e(f, current_time_step);
        _mz_xn(f, 0, j - js, k - ks) +=
            mz_xn * _transform_e(f, current_time_step);
        _jy_xn(f, 0, j - js, k - ks) +=
            jy_xn * _transform_h(f, current_time_step);
        _jz_xn(f, 0, j - js, k - ks) +=
            jz_xn * _transform_h(f, current_time_step);
        _my_xp(f, 0, j - js, k - ks) +=
            my_xp * _transform_e(f, current_time_step);
        _mz_xp(f, 0, j - js, k - ks) +=
            mz_xp * _transform_e(f, current_time_step);
        _jy_xp(f, 0, j - js, k - ks) +=
            jy_xp * _transform_h(f, current_time_step);
        _jz_xp(f, 0, j - js, k - ks) +=
            jz_xp * _transform_h(f, current_time_step);
      }
    }
  }

  for (auto i{is}; i < ie; ++i) {
    for (auto k{ks}; k < ke; ++k) {
      auto mz_yn{-1.0 * emf->exFaceCenterY(i, js, k)};
      auto mx_yn{emf->ezFaceCenterY(i, js, k)};
      auto jz_yn{emf->hxFaceCenterY(i, js, k)};
      auto jx_yn{-1.0 * emf->hzFaceCenterY(i, js, k)};
      auto mz_yp{emf->exFaceCenterY(i, je, k)};
      auto mx_yp{-1.0 * emf->ezFaceCenterY(i, je, k)};
      auto jz_yp{-1.0 * emf->hxFaceCenterY(i, je, k)};
      auto jx_yp{emf->hzFaceCenterY(i, je, k)};

      for (int f{0}; f < _frequencies.size(); ++f) {
        _mx_yn(f, i - is, 0, k - ks) +=
            mx_yn * _transform_e(f, current_time_step);
        _mz_yn(f, i - is, 0, k - ks) +=
            mz_yn * _transform_e(f, current_time_step);
        _jx_yn(f, i - is, 0, k - ks) +=
            jx_yn * _transform_h(f, current_time_step);
        _jz_yn(f, i - is, 0, k - ks) +=
            jz_yn * _transform_h(f, current_time_step);
        _mx_yp(f, i - is, 0, k - ks) +=
            mx_yp * _transform_e(f, current_time_step);
        _mz_yp(f, i - is, 0, k - ks) +=
            mz_yp * _transform_e(f, current_time_step);
        _jx_yp(f, i - is, 0, k - ks) +=
            jx_yp * _transform_h(f, current_time_step);
        _jz_yp(f, i - is, 0, k - ks) +=
            jz_yp * _transform_h(f, current_time_step);
      }
    }
  }

  for (auto i{is}; i < ie; ++i) {
    for (auto j{js}; j < js; ++j) {
      auto mx_zn{-1.0 * emf->eyFaceCenterZ(i, j, ks)};
      auto my_zn{emf->exFaceCenterZ(i, j, ks)};
      auto jx_zn{emf->hyFaceCenterZ(i, j, ks)};
      auto jy_zn{-1.0 * emf->hxFaceCenterZ(i, j, ks)};
      auto mx_zp{emf->eyFaceCenterZ(i, j, ke)};
      auto my_zp{-1.0 * emf->exFaceCenterZ(i, j, ke)};
      auto jx_zp{-1.0 * emf->hyFaceCenterZ(i, j, ke)};
      auto jy_zp{emf->hxFaceCenterZ(i, j, ke)};
      for (int f{0}; f < _frequencies.size(); ++f) {
        _mx_zn(f, i - is, j - js, 0) +=
            mx_zn * _transform_e(f, current_time_step);
        _my_zn(f, i - is, j - js, 0) +=
            my_zn * _transform_e(f, current_time_step);
        _jx_zn(f, i - is, j - js, 0) +=
            jx_zn * _transform_h(f, current_time_step);
        _jy_zn(f, i - is, j - js, 0) +=
            jy_zn * _transform_h(f, current_time_step);
        _mx_zp(f, i - is, j - js, 0) +=
            mx_zp * _transform_e(f, current_time_step);
        _my_zp(f, i - is, j - js, 0) +=
            my_zp * _transform_e(f, current_time_step);
        _jx_zp(f, i - is, j - js, 0) +=
            jx_zp * _transform_h(f, current_time_step);
        _jy_zp(f, i - is, j - js, 0) +=
            jy_zp * _transform_h(f, current_time_step);
      }
    }
  }
}

void NFFFT::output() {}

void NFFFT::processFarField(const xt::xarray<double>& theta,
                            const xt::xarray<double>& phi,
                            const Vector& origin = Vector{0.0, 0.0, 0.0}) {
  using namespace std::complex_literals;
  auto box{gridBoxPtr()};
  const auto is{box->origin().i()};
  const auto js{box->origin().j()};
  const auto ks{box->origin().k()};
  const auto ie{box->end().i()};
  const auto je{box->end().j()};
  const auto ke{box->end().k()};

  auto g{gridSpacePtr()};
  const auto& size_x{_grid_space->eSizeX()};
  const auto& size_y{_grid_space->eSizeY()};
  const auto& size_z{_grid_space->eSizeZ()};

  const auto& h_node_x{g->hNodeX()};
  const auto& h_node_y{g->hNodeY()};
  const auto& h_node_z{g->hNodeZ()};
  const auto& e_node_x{g->eNodeX()};
  const auto& e_node_y{g->eNodeY()};
  const auto& e_node_z{g->eNodeZ()};

  _a_theta = xt::zeros<std::complex<double>>(
      {_frequencies.size(), theta.size(), phi.size()});
  _a_phi = xt::zeros<std::complex<double>>(
      {_frequencies.size(), theta.size(), phi.size()});
  _f_theta = xt::zeros<std::complex<double>>(
      {_frequencies.size(), theta.size(), phi.size()});
  _f_phi = xt::zeros<std::complex<double>>(
      {_frequencies.size(), theta.size(), phi.size()});

  for (auto f{0}; f < _frequencies.size(); ++f) {
    auto wave_number{2 * constant::PI * _frequencies(f) / constant::C_0};
    _a_x = 0.0 + 0.0i;
    _a_y = 0.0 + 0.0i;
    _a_z = 0.0 + 0.0i;
    _f_x = 0.0 + 0.0i;
    _f_y = 0.0 + 0.0i;
    _f_z = 0.0 + 0.0i;

    for (std::size_t t{0}; t < theta.size(); ++t) {
      for (std::size_t p{0}; p < phi.size(); ++p) {
        auto sin_t_cos_p{std::sin(theta(t)) * std::cos(phi(p))};
        auto sin_t_sin_p{std::sin(theta(t)) * std::sin(phi(p))};
        auto cos_t{std::cos(theta(t))};

        for (auto j{js}; j < je; ++j) {
          for (auto k{ks}; k < ke; ++k) {
            auto ds{size_y(j) * size_z(k)};
            auto delay_xn{std::abs(origin.x() - e_node_x(is)) * sin_t_cos_p +
                          std::abs(origin.y() - h_node_y(j)) * sin_t_sin_p +
                          std::abs(origin.z() - h_node_z(k)) * cos_t};
            auto phase_shift_xn{std::exp(1i * wave_number * delay_xn)};
            auto delay_xp{std::abs(origin.x() - e_node_x(ie)) * sin_t_cos_p +
                          std::abs(origin.y() - h_node_y(j)) * sin_t_sin_p +
                          std::abs(origin.z() - h_node_z(k)) * cos_t};
            auto phase_shift_xp{std::exp(1i * wave_number * delay_xp)};
            _a_y += ds * phase_shift_xn * _jy_xn(f, 0, j - js, k - ks);
            _a_z += ds * phase_shift_xn * _jz_xn(f, 0, j - js, k - ks);
            _f_y += ds * phase_shift_xn * _my_xn(f, 0, j - js, k - ks);
            _f_z += ds * phase_shift_xn * _mz_xn(f, 0, j - js, k - ks);
            _a_y += ds * phase_shift_xp * _jy_xp(f, 0, j - js, k - ks);
            _a_z += ds * phase_shift_xp * _jz_xp(f, 0, j - js, k - ks);
            _f_y += ds * phase_shift_xp * _my_xp(f, 0, j - js, k - ks);
            _f_z += ds * phase_shift_xp * _mz_xp(f, 0, j - js, k - ks);
          }
        }

        for (auto i{is}; i < ie; ++i) {
          for (auto k{ks}; k < ke; ++k) {
            auto ds{size_x(i) * size_z(k)};
            auto delay_yn{std::abs(origin.x() - h_node_x(i)) * sin_t_cos_p +
                          std::abs(origin.y() - e_node_y(js)) * sin_t_sin_p +
                          std::abs(origin.z() - h_node_z(k)) * cos_t};
            auto phase_shift_yn{std::exp(1i * wave_number * delay_yn)};
            auto delay_yp{std::abs(origin.x() - h_node_x(i)) * sin_t_cos_p +
                          std::abs(origin.y() - e_node_y(je)) * sin_t_sin_p +
                          std::abs(origin.z() - h_node_z(k)) * cos_t};
            auto phase_shift_yp{std::exp(1i * wave_number * delay_yp)};
            _a_x += ds * phase_shift_yn * _jx_yn(f, i - is, 0, k - ks);
            _a_z += ds * phase_shift_yn * _jz_yn(f, i - is, 0, k - ks);
            _f_x += ds * phase_shift_yn * _mx_yn(f, i - is, 0, k - ks);
            _f_z += ds * phase_shift_yn * _mz_yn(f, i - is, 0, k - ks);
            _a_x += ds * phase_shift_yp * _jx_yp(f, i - is, 0, k - ks);
            _a_z += ds * phase_shift_yp * _jz_yp(f, i - is, 0, k - ks);
            _f_x += ds * phase_shift_yp * _mx_yp(f, i - is, 0, k - ks);
            _f_z += ds * phase_shift_yp * _mz_yp(f, i - is, 0, k - ks);
          }
        }

        for (auto i{is}; i < ie; ++i) {
          for (auto j{js}; j < je; ++j) {
            auto ds{size_x(i) * size_y(j)};
            auto delay_zn{std::abs(origin.x() - h_node_x(i)) * sin_t_cos_p +
                          std::abs(origin.y() - h_node_y(j)) * sin_t_sin_p +
                          std::abs(origin.z() - e_node_z(ks)) * cos_t};
            auto phase_shift_zn{std::exp(1i * wave_number * delay_zn)};
            auto delay_zp{std::abs(origin.x() - h_node_x(i)) * sin_t_cos_p +
                          std::abs(origin.y() - h_node_y(j)) * sin_t_sin_p +
                          std::abs(origin.z() - e_node_z(ke)) * cos_t};
            auto phase_shift_zp{std::exp(1i * wave_number * delay_zp)};
            _a_x += ds * phase_shift_zn * _jx_zn(f, i - is, j - js, 0);
            _a_y += ds * phase_shift_zn * _jy_zn(f, i - is, j - js, 0);
            _f_x += ds * phase_shift_zn * _mx_zn(f, i - is, j - js, 0);
            _f_y += ds * phase_shift_zn * _my_zn(f, i - is, j - js, 0);
            _a_x += ds * phase_shift_zp * _jx_zp(f, i - is, j - js, 0);
            _a_y += ds * phase_shift_zp * _jy_zp(f, i - is, j - js, 0);
            _f_x += ds * phase_shift_zp * _mx_zp(f, i - is, j - js, 0);
            _f_y += ds * phase_shift_zp * _my_zp(f, i - is, j - js, 0);
          }
        }

        auto sin_t{std::sin(theta(t))};
        auto cos_p{std::cos(p)};
        auto sin_p{std::sin(p)};
        auto cos_t_cos_p{cos_t * cos_p};
        auto cos_t_sin_p{cos_t * sin_p};

        _a_theta(f, t, p) =
            cos_t_cos_p * _a_x + cos_t_sin_p * _a_y - sin_t * _a_z;
        _a_phi(f, t, p) = -sin_p * _a_x + cos_p * _a_y;
        _f_theta(f, t, p) =
            cos_t_cos_p * _f_x + cos_t_sin_p * _f_y - sin_t * _f_z;
        _f_phi(f, t, p) = -sin_p * _f_x + cos_p * _f_y;
      }
    }
  }

  auto out_dir{std::filesystem::path{_output_dir}};
  if (!std::filesystem::exists(out_dir) ||
      !std::filesystem::is_directory(out_dir)) {
    std::filesystem::create_directories(out_dir);
  }

  auto a_theta_path{out_dir / "a_theta.npy"};
  auto a_phi_path{out_dir / "a_phi.npy"};
  auto f_theta_path{out_dir / "f_theta.npy"};
  auto f_phi_path{out_dir / "f_phi.npy"};

  xt::dump_npy(a_theta_path.string(), _a_theta);
  xt::dump_npy(a_phi_path.string(), _a_phi);
  xt::dump_npy(f_theta_path.string(), _f_theta);
  xt::dump_npy(f_phi_path.string(), _f_phi);
}

}  // namespace xfdtd

#endif  // _XFDTD_LIB_NFFFT_H_
