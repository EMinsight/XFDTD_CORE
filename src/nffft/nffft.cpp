#include <xfdtd/nffft/nffft.h>
#include <xfdtd/util/constant.h>

#include <complex>
#include <cstdlib>
#include <utility>
#include <xtensor.hpp>
#include <xtensor/xnpy.hpp>

namespace xfdtd {

NFFFT::NFFFT(std::size_t distance_x, std::size_t distance_y,
             std::size_t distance_z, xt::xarray<double> frequencies,
             std::string output_dir)
    : _distance_x{distance_x},
      _distance_y{distance_y},
      _distance_z{distance_z},
      _frequencies{std::move(frequencies)},
      _output_dir{std::move(output_dir)} {}

const GridSpace* NFFFT::gridSpacePtr() const { return _grid_space.get(); }

const CalculationParam* NFFFT::calculationParamPtr() const {
  return _calculation_param.get();
}

const EMF* NFFFT::emfPtr() const { return _emf.get(); }

const GridBox* NFFFT::gridBoxPtr() const { return _grid_box.get(); }

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
      _transform_e(f, t) = dt * std::exp(-1i * 2.0 * constant::PI *
                                         (_frequencies(f) * (t + 1) * dt));
      _transform_h(f, t) = dt * std::exp(-1i * 2.0 * constant::PI *
                                         (_frequencies(f) * (t + 0.5) * dt));
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

      auto my_xp{0.5 * (emf->ez()(ie, j + 1, k) + emf->ez()(ie, j, k))};
      auto mz_xp{-0.5 * (emf->ey()(ie, j, k + 1) + emf->ey()(ie, j, k))};
      auto jy_xp{-0.25 *
                 (emf->hz()(ie, j, k + 1) + emf->hz()(ie, j, k) +
                  emf->hz()(ie - 1, j, k + 1) + emf->hz()(ie - 1, j, k))};
      auto jz_xp{0.25 *
                 (emf->hy()(ie, j + 1, k) + emf->hy()(ie, j, k) +
                  emf->hy()(ie - 1, j + 1, k) + emf->hy()(ie - 1, j, k))};

      for (std::size_t f{0}; f < _frequencies.size(); ++f) {
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
      auto mz_yn{-0.5 * (emf->ex()(i, js, k + 1) + emf->ex()(i, js, k))};
      auto mx_yn{0.5 * (emf->ez()(i + 1, js, k) + emf->ez()(i, js, k))};
      auto jz_yn{0.25 *
                 (emf->hx()(i + 1, js, k) + emf->hx()(i, js, k) +
                  emf->hx()(i + 1, js - 1, k) + emf->hx()(i, js - 1, k))};
      auto jx_yn{-0.25 *
                 (emf->hz()(i, js, k + 1) + emf->hz()(i, js, k) +
                  emf->hz()(i, js - 1, k + 1) + emf->hz()(i, js - 1, k))};

      auto mz_yp{0.5 * (emf->ex()(i, je, k + 1) + emf->ex()(i, je, k))};
      auto mx_yp{-0.5 * (emf->ez()(i + 1, je, k) + emf->ez()(i, je, k))};
      auto jz_yp{-0.25 *
                 (emf->hx()(i + 1, je, k) + emf->hx()(i, je, k) +
                  emf->hx()(i + 1, je - 1, k) + emf->hx()(i, je - 1, k))};
      auto jx_yp{0.25 *
                 (emf->hz()(i, je, k + 1) + emf->hz()(i, je, k) +
                  emf->hz()(i, je - 1, k + 1) + emf->hz()(i, je - 1, k))};

      for (std::size_t f{0}; f < _frequencies.size(); ++f) {
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
    for (auto j{js}; j < je; ++j) {
      auto mx_zn{-0.5 * (emf->ey()(i + 1, j, ks) + emf->ey()(i, j, ks))};
      auto my_zn{0.5 * (emf->ex()(i, j + 1, ks) + emf->ex()(i, j, ks))};
      auto jx_zn{0.25 *
                 (emf->hy()(i, j + 1, ks) + emf->hy()(i, j, ks) +
                  emf->hy()(i, j + 1, ks - 1) + emf->hy()(i, j, ks - 1))};
      auto jy_zn{-0.25 *
                 (emf->hx()(i + 1, j, ks) + emf->hx()(i, j, ks) +
                  emf->hx()(i + 1, j, ks - 1) + emf->hx()(i, j, ks - 1))};

      auto mx_zp{0.5 * (emf->ey()(i + 1, j, ke) + emf->ey()(i, j, ke))};
      auto my_zp{-0.5 * (emf->ex()(i, j + 1, ke) + emf->ex()(i, j, ke))};
      auto jx_zp{-0.25 *
                 (emf->hy()(i, j + 1, ke) + emf->hy()(i, j, ke) +
                  emf->hy()(i, j + 1, ke - 1) + emf->hy()(i, j, ke - 1))};
      auto jy_zp{0.25 *
                 (emf->hx()(i + 1, j, ke) + emf->hx()(i, j, ke) +
                  emf->hx()(i + 1, j, ke - 1) + emf->hx()(i, j, ke - 1))};

      for (std::size_t f{0}; f < _frequencies.size(); ++f) {
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
                            const std::string& sub_dir, const Vector& origin) {
  using namespace std::complex_literals;
  auto box{gridBoxPtr()};
  const auto is{box->origin().i()};
  const auto js{box->origin().j()};
  const auto ks{box->origin().k()};
  const auto ie{box->end().i()};
  const auto je{box->end().j()};
  const auto ke{box->end().k()};

  auto g{gridSpacePtr()};
  const auto& size_x{g->eSizeX()};
  const auto& size_y{g->eSizeY()};
  const auto& size_z{g->eSizeZ()};

  const auto& h_node_x{g->hNodeX()};
  const auto& h_node_y{g->hNodeY()};
  const auto& h_node_z{g->hNodeZ()};
  const auto& e_node_x{g->eNodeX()};
  const auto& e_node_y{g->eNodeY()};
  const auto& e_node_z{g->eNodeZ()};

  if (theta.size() != 1 && phi.size() != 1) {
    throw std::runtime_error("theta.size() != 1 || phi.size() != 0");
  }

  auto num_angle{theta.size() * phi.size()};

  _a_theta = xt::zeros<std::complex<double>>({_frequencies.size(), num_angle});
  _a_phi = xt::zeros<std::complex<double>>({_frequencies.size(), num_angle});
  _f_theta = xt::zeros<std::complex<double>>({_frequencies.size(), num_angle});
  _f_phi = xt::zeros<std::complex<double>>({_frequencies.size(), num_angle});

  auto cos_t{xt::cos(theta)};
  auto sin_t{xt::sin(theta)};
  auto cos_p{xt::cos(phi)};
  auto sin_p{xt::sin(phi)};
  auto sin_t_sin_p{sin_t * sin_p};
  auto sin_t_cos_p{sin_t * cos_p};
  auto cos_t_sin_p{cos_t * sin_p};
  auto cos_t_cos_p{cos_t * cos_p};

  for (std::size_t f{0}; f < _frequencies.size(); ++f) {
    auto wave_number{2.0 * constant::PI * _frequencies(f) / constant::C_0};

    for (auto j{js}; j < je; ++j) {
      for (auto k{ks}; k < ke; ++k) {
        auto ds{size_y(j) * size_z(k)};
        auto r_xn{Vector{e_node_x(is), h_node_y(j), h_node_z(k)} - origin};
        auto phase_shift_xn{
            xt::exp(1i * wave_number *
                    (r_xn.x() * sin_t_cos_p + r_xn.y() * sin_t_sin_p +
                     r_xn.z() * cos_t))};
        xt::view(_a_theta, f, xt::all()) +=
            (_jy_xn(f, 0, j - js, k - ks) * cos_t_sin_p -
             _jz_xn(f, 0, j - js, k - ks) * sin_t) *
            phase_shift_xn * ds;
        xt::view(_f_theta, f, xt::all()) +=
            (_my_xn(f, 0, j - js, k - ks) * cos_t_sin_p -
             _mz_xn(f, 0, j - js, k - ks) * sin_t) *
            phase_shift_xn * ds;
        xt::view(_a_phi, f, xt::all()) +=
            (_jy_xn(f, 0, j - js, k - ks) * cos_p) * phase_shift_xn * ds;
        xt::view(_f_phi, f, xt::all()) +=
            (_my_xn(f, 0, j - js, k - ks) * cos_p) * phase_shift_xn * ds;
        auto r_xp{Vector{e_node_x(ie), h_node_y(j), h_node_z(k)} - origin};
        auto phase_shift_xp{
            xt::exp(1i * wave_number *
                    (r_xp.x() * sin_t_cos_p + r_xp.y() * sin_t_sin_p +
                     r_xp.z() * cos_t))};
        xt::view(_a_theta, f, xt::all()) +=
            (_jy_xp(f, 0, j - js, k - ks) * cos_t_sin_p -
             _jz_xp(f, 0, j - js, k - ks) * sin_t) *
            phase_shift_xp * ds;
        xt::view(_f_theta, f, xt::all()) +=
            (_my_xp(f, 0, j - js, k - ks) * cos_t_sin_p -
             _mz_xp(f, 0, j - js, k - ks) * sin_t) *
            phase_shift_xp * ds;
        xt::view(_a_phi, f, xt::all()) +=
            (_jy_xp(f, 0, j - js, k - ks) * cos_p) * phase_shift_xp * ds;
        xt::view(_f_phi, f, xt::all()) +=
            (_my_xp(f, 0, j - js, k - ks) * cos_p) * phase_shift_xp * ds;
      }
    }

    for (auto i{is}; i < ie; ++i) {
      for (auto k{ks}; k < ke; ++k) {
        auto ds{size_x(i) * size_z(k)};
        auto r_yn{Vector{h_node_x(i), e_node_y(js), h_node_z(k)} - origin};
        auto phase_shift_yn{
            xt::exp(1i * wave_number *
                    (r_yn.x() * sin_t_cos_p + r_yn.y() * sin_t_sin_p +
                     r_yn.z() * cos_t))};
        xt::view(_a_theta, f, xt::all()) +=
            (_jx_yn(f, i - is, 0, k - ks) * cos_t_cos_p -
             _jz_yn(f, i - is, 0, k - ks) * sin_t) *
            phase_shift_yn * ds;
        xt::view(_f_theta, f, xt::all()) +=
            (_mx_yn(f, i - is, 0, k - ks) * cos_t_cos_p -
             _mz_yn(f, i - is, 0, k - ks) * sin_t) *
            phase_shift_yn * ds;
        xt::view(_a_phi, f, xt::all()) -=
            (_jx_yn(f, i - is, 0, k - ks) * sin_p) * phase_shift_yn * ds;
        xt::view(_f_phi, f, xt::all()) -=
            (_mx_yn(f, i - is, 0, k - ks) * sin_p) * phase_shift_yn * ds;
        auto r_yp{Vector{h_node_x(i), e_node_y(je), h_node_z(k)} - origin};
        auto phase_shift_yp{
            xt::exp(1i * wave_number *
                    (r_yp.x() * sin_t_cos_p + r_yp.y() * sin_t_sin_p +
                     r_yp.z() * cos_t))};
        xt::view(_a_theta, f, xt::all()) +=
            (_jx_yp(f, i - is, 0, k - ks) * cos_t_cos_p -
             _jz_yp(f, i - is, 0, k - ks) * sin_t) *
            phase_shift_yp * ds;
        xt::view(_f_theta, f, xt::all()) +=
            (_mx_yp(f, i - is, 0, k - ks) * cos_t_cos_p -
             _mz_yp(f, i - is, 0, k - ks) * sin_t) *
            phase_shift_yp * ds;
        xt::view(_a_phi, f, xt::all()) -=
            (_jx_yp(f, i - is, 0, k - ks) * sin_p) * phase_shift_yp * ds;
        xt::view(_f_phi, f, xt::all()) -=
            (_mx_yp(f, i - is, 0, k - ks) * sin_p) * phase_shift_yp * ds;
      }
    }

    for (auto i{is}; i < ie; ++i) {
      for (auto j{js}; j < je; ++j) {
        auto ds = size_x(i) * size_y(j);
        auto r_zn{Vector{h_node_x(i), h_node_y(j), e_node_z(ks)} - origin};
        auto phase_shift_zn{
            xt::exp(1i * wave_number *
                    (r_zn.x() * sin_t_cos_p + r_zn.y() * sin_t_sin_p +
                     r_zn.z() * cos_t))};
        xt::view(_a_theta, f, xt::all()) +=
            (_jx_zn(f, i - is, j - js, 0) * cos_t_cos_p +
             _jy_zn(f, i - is, j - js, 0) * cos_t_sin_p) *
            phase_shift_zn * ds;
        xt::view(_f_theta, f, xt::all()) +=
            (_mx_zn(f, i - is, j - js, 0) * cos_t_cos_p +
             _my_zn(f, i - is, j - js, 0) * cos_t_sin_p) *
            phase_shift_zn * ds;
        xt::view(_a_phi, f, xt::all()) +=
            (-_jx_zn(f, i - is, j - js, 0) * sin_p +
             _jy_zn(f, i - is, j - js, 0) * cos_p) *
            phase_shift_zn * ds;
        xt::view(_f_phi, f, xt::all()) +=
            (-_mx_zn(f, i - is, j - js, 0) * sin_p +
             _my_zn(f, i - is, j - js, 0) * cos_p) *
            phase_shift_zn * ds;
        auto r_zp{Vector{h_node_x(i), h_node_y(j), e_node_z(ke)} - origin};
        auto phase_shift_zp{
            xt::exp(1i * wave_number *
                    (r_zp.x() * sin_t_cos_p + r_zp.y() * sin_t_sin_p +
                     r_zp.z() * cos_t))};
        xt::view(_a_theta, f, xt::all()) +=
            (_jx_zp(f, i - is, j - js, 0) * cos_t_cos_p +
             _jy_zp(f, i - is, j - js, 0) * cos_t_sin_p) *
            phase_shift_zp * ds;
        xt::view(_f_theta, f, xt::all()) +=
            (_mx_zp(f, i - is, j - js, 0) * cos_t_cos_p +
             _my_zp(f, i - is, j - js, 0) * cos_t_sin_p) *
            phase_shift_zp * ds;
        xt::view(_a_phi, f, xt::all()) +=
            (-_jx_zp(f, i - is, j - js, 0) * sin_p +
             _jy_zp(f, i - is, j - js, 0) * cos_p) *
            phase_shift_zp * ds;
        xt::view(_f_phi, f, xt::all()) +=
            (-_mx_zp(f, i - is, j - js, 0) * sin_p +
             _my_zp(f, i - is, j - js, 0) * cos_p) *
            phase_shift_zp * ds;
      }
    }
  }

  auto out_dir{std::filesystem::path{_output_dir} / sub_dir};
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

void NFFFT::outputRadiationPower() {
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

  xt::xarray<double> res{xt::zeros<double>({_frequencies.size()})};

  for (std::size_t f{0}; f < _frequencies.size(); ++f) {
    auto power{std::complex<double>{0.0, 0.0}};
    auto dy_dz{size_y(0) * size_z(0)};
    auto dz_dx{size_z(0) * size_x(0)};
    auto dx_dy{size_x(0) * size_y(0)};

    for (std::size_t j{js}; j < je; ++j) {
      for (std::size_t k{ks}; k < ke; ++k) {
        auto ds{size_y(j) * size_z(k)};
        power -= ds * (_mz_xn(f, 0, j - js, k - ks) *
                           std::conj(_jy_xn(f, 0, j - js, k - ks)) -
                       _my_xn(f, 0, j - js, k - ks) *
                           std::conj(_jz_xn(f, 0, j - js, k - ks)));
        power += ds * (_mz_xp(f, 0, j - js, k - ks) *
                           std::conj(_jy_xp(f, 0, j - js, k - ks)) -
                       _my_xp(f, 0, j - js, k - ks) *
                           std::conj(_jz_xp(f, 0, j - js, k - ks)));
      }
    }
    for (std::size_t i{is}; i < ie; ++i) {
      for (std::size_t k{ks}; k < ke; ++k) {
        auto ds{size_x(i) * size_z(k)};
        power -= ds * (_mx_yn(f, i - is, 0, k - ks) *
                           std::conj(_jx_yn(f, i - is, 0, k - ks)) -
                       _mz_yn(f, i - is, 0, k - ks) *
                           std::conj(_jz_yn(f, i - is, 0, k - ks)));
        power += ds * (_mx_yp(f, i - is, 0, k - ks) *
                           std::conj(_jx_yp(f, i - is, 0, k - ks)) -
                       _mz_yp(f, i - is, 0, k - ks) *
                           std::conj(_jz_yp(f, i - is, 0, k - ks)));
      }
    }

    for (std::size_t i{is}; i < ie; ++i) {
      for (std::size_t j{js}; j < je; ++j) {
        auto ds{size_x(i) * size_y(j)};
        power -= ds * (_my_zn(f, i - is, j - js, 0) *
                           std::conj(_jy_zn(f, i - is, j - js, 0)) -
                       _mx_zn(f, i - is, j - js, 0) *
                           std::conj(_jx_zn(f, i - is, j - js, 0)));
        power += ds * (_my_zp(f, i - is, j - js, 0) *
                           std::conj(_jy_zp(f, i - is, j - js, 0)) -
                       _mx_zp(f, i - is, j - js, 0) *
                           std::conj(_jx_zp(f, i - is, j - js, 0)));
      }
    }

    res(f) = 0.5 * std::real(power);
  }

  auto radiation_power_path{std::filesystem::path{_output_dir} /
                            "radiation_power.npy"};
  if (!std::filesystem::exists(radiation_power_path) ||
      !std::filesystem::is_regular_file(radiation_power_path)) {
    std::filesystem::create_directories(radiation_power_path.parent_path());
  }
  xt::dump_npy(radiation_power_path.string(), res);
}

}  // namespace xfdtd
