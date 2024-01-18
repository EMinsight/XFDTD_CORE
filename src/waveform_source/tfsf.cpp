#include <xfdtd/waveform_source/tfsf.h>

#include <cmath>
#include <cstdlib>
#include <xtensor-blas/xlinalg.hpp>

#include "xfdtd/util/constant.h"

namespace xfdtd {

TFSF::TFSF(std::size_t x, std::size_t y, std::size_t z, double theta,
           double phi, double psi, std::unique_ptr<Waveform> waveform)
    : WaveformSource{std::move(waveform)},
      _x{x},
      _y{y},
      _z{z},
      _theta{theta},
      _phi{phi},
      _psi{psi},
      _sin_theta{std::sin(theta)},
      _cos_theta{std::cos(theta)},
      _sin_phi{std::sin(phi)},
      _cos_phi{std::cos(phi)},
      _sin_psi{std::sin(psi)},
      _cos_psi{std::cos(psi)} {
  initTransform();
}

void TFSF::correctMaterialSpace() {}

void TFSF::correctUpdateCoefficient() {}

void TFSF::updateWaveformSourceE() {
  for (std::size_t i{0}; i < _h_inc.size(); ++i) {
    _h_inc(i) = _chih * _h_inc(i) + _chiei * (_e_inc(i + 1) - _e_inc(i));
  }
  auto x{_e_inc(_e_inc.size() - 2)};
  auto y{_e_inc(_e_inc.size() - 1)};
  _e_inc(0) = waveform()->value()(
      calculationParamPtr()->timeParam()->currentTimeStep());
  for (std::size_t i{1}; i < _e_inc.size(); ++i) {
    _e_inc(i) = _ceie * _e_inc(i) + _ceihi * (_h_inc(i) - _h_inc(i - 1));
  }

  _e_inc(_e_inc.size() - 1) = -_a +
                              _abc_coff_0 * (_e_inc(_e_inc.size() - 2) + _b) +
                              _abc_coff_1 * (x + y);
  _a = x;
  _b = y;

  _ex_inc = _transform_e(0) * _e_inc;
  _ey_inc = _transform_e(1) * _e_inc;
  _ez_inc = _transform_e(2) * _e_inc;
  _hx_inc = _transform_h(0) * _h_inc;
  _hy_inc = _transform_h(1) * _h_inc;
  _hz_inc = _transform_h(2) * _h_inc;
}

void TFSF::updateWaveformSourceH() {}

std::size_t TFSF::x() const { return _x; }

std::size_t TFSF::y() const { return _y; }

std::size_t TFSF::z() const { return _z; }

double TFSF::theta() const { return _theta; }

double TFSF::phi() const { return _phi; }

double TFSF::psi() const { return _psi; }

double TFSF::sinTheta() const { return _sin_theta; }

double TFSF::cosTheta() const { return _cos_theta; }

double TFSF::sinPhi() const { return _sin_phi; }

double TFSF::cosPhi() const { return _cos_phi; }

double TFSF::sinPsi() const { return _sin_psi; }

double TFSF::cosPsi() const { return _cos_psi; }

Vector TFSF::k() const { return _k; }

const GridBox *TFSF::boxPtr() const { return _box.get(); }

void TFSF::defaultInit(std::shared_ptr<GridSpace> grid_space,
                       std::shared_ptr<CalculationParam> calculation_param,
                       std::shared_ptr<EMF> emf) {
  WaveformSource::defaultInit(std::move(grid_space),
                              std::move(calculation_param), std::move(emf));

  if (gridSpacePtr()->type() != GridSpace::Type::UNIFORM) {
    throw std::runtime_error("TFSF only supports uniform grid space");
  }

  auto origin_x{gridSpacePtr()->box().origin().i() + x()};
  auto origin_y{gridSpacePtr()->box().origin().j() + y()};
  auto origin_z{gridSpacePtr()->box().origin().k() + z()};
  auto size_x{gridSpacePtr()->box().size().i() - 2 * x()};
  auto size_y{gridSpacePtr()->box().size().j() - 2 * y()};
  auto size_z{gridSpacePtr()->box().size().k() - 2 * z()};
  _box = std::make_unique<GridBox>(Grid{origin_x, origin_y, origin_z},
                                   Grid{size_x, size_y, size_z});

  waveform()->init(calculationParamPtr()->timeParam()->eTime() -
                   calculationParamPtr()->timeParam()->dt());
  _ratio_delta =
      1 / std::sqrt(std::pow(sinTheta(), 4) *
                        (std::pow(cosPhi(), 4) + std::pow(sinPhi(), 4)) +
                    std::pow(cosTheta(), 4));
  _auxiliary_size =
      std::ceil<std::size_t>(
          _ratio_delta *
          (std::sqrt(pow(size_x, 2) + pow(size_y, 2) + pow(size_z, 2)))) +
      4 + 1;

  calculateProjection();

  _e_inc = xt::zeros<double>({_auxiliary_size});
  _ex_inc = xt::zeros<double>({_auxiliary_size});
  _ey_inc = xt::zeros<double>({_auxiliary_size});
  _ez_inc = xt::zeros<double>({_auxiliary_size});
  _h_inc = xt::zeros<double>({_auxiliary_size - 1});
  _hx_inc = xt::zeros<double>({_auxiliary_size - 1});
  _hy_inc = xt::zeros<double>({_auxiliary_size - 1});
  _hz_inc = xt::zeros<double>({_auxiliary_size - 1});

  if (gridSpacePtr()->dimension() == GridSpace::Dimension::ONE) {
    _scaled_dl = gridSpacePtr()->basedDz() / _ratio_delta;
  }
  if (gridSpacePtr()->dimension() == GridSpace::Dimension::TWO) {
    _scaled_dl = gridSpacePtr()->basedDx() / _ratio_delta;
  }
  if (gridSpacePtr()->dimension() == GridSpace::Dimension::THREE) {
    _scaled_dl = gridSpacePtr()->basedDz() / _ratio_delta;
  }

  _ceie = 1;
  _chih = 1;
  _ceihi = -(calculationParamPtr()->timeParam()->dt() /
             (constant::EPSILON_0 * _scaled_dl));
  _chiei = -(calculationParamPtr()->timeParam()->dt() /
             (constant::MU_0 * _scaled_dl));
  _abc_coff_0 =
      (constant::C_0 * calculationParamPtr()->timeParam()->dt() - _scaled_dl) /
      (constant::C_0 * calculationParamPtr()->timeParam()->dt() + _scaled_dl);
  _abc_coff_1 =
      2 * _scaled_dl /
      (constant::C_0 * calculationParamPtr()->timeParam()->dt() + _scaled_dl);
}

double TFSF::exInc(std::size_t i, std::size_t j, std::size_t k) {
  i = i - boxPtr()->origin().i() + 1;
  j = j - boxPtr()->origin().j();
  k = k - boxPtr()->origin().k();
  auto projection{_projection_x_half(i) + _projection_y_int(j) +
                  _projection_z_int(k)};
  auto index{static_cast<std::size_t>(projection)};
  auto weight{projection - index};
  return (1 - weight) * _ex_inc(index) + weight * _ex_inc(index + 1);
}

double TFSF::eyInc(std::size_t i, std::size_t j, std::size_t k) {
  i = i - boxPtr()->origin().i();
  j = j - boxPtr()->origin().j() + 1;
  k = k - boxPtr()->origin().k();
  auto projection{_projection_x_int(i) + _projection_y_half(j) +
                  _projection_z_int(k)};
  auto index{static_cast<std::size_t>(projection)};
  auto weight{projection - index};
  return (1 - weight) * _ey_inc(index) + weight * _ey_inc(index + 1);
}

double TFSF::ezInc(std::size_t i, std::size_t j, std::size_t k) {
  i = i - boxPtr()->origin().i();
  j = j - boxPtr()->origin().j();
  k = k - boxPtr()->origin().k() + 1;
  auto projection{_projection_x_int(i) + _projection_y_int(j) +
                  _projection_z_half(k)};
  auto index{static_cast<std::size_t>(projection)};
  auto weight{projection - index};
  return (1 - weight) * _ez_inc(index) + weight * _ez_inc(index + 1);
}

double TFSF::hxInc(std::size_t i, std::size_t j, std::size_t k) {
  i = i - boxPtr()->origin().i();
  j = j - boxPtr()->origin().j() + 1;
  k = k - boxPtr()->origin().k() + 1;
  auto projection{_projection_x_int(i) + _projection_y_half(j) +
                  _projection_z_half(k) - 0.5};
  auto index{static_cast<std::size_t>(projection)};
  auto weight{projection - index};
  return (1 - weight) * _hx_inc(index) + weight * _hx_inc(index + 1);
}

double TFSF::hyInc(std::size_t i, std::size_t j, std::size_t k) {
  i = i - boxPtr()->origin().i() + 1;
  j = j - boxPtr()->origin().j();
  k = k - boxPtr()->origin().k() + 1;
  auto projection{_projection_x_half(i) + _projection_y_int(j) +
                  _projection_z_half(k) - 0.5};
  auto index{static_cast<std::size_t>(projection)};
  auto weight{projection - index};
  return (1 - weight) * _hy_inc(index) + weight * _hy_inc(index + 1);
}

double TFSF::hzInc(std::size_t i, std::size_t j, std::size_t k) {
  i = i - boxPtr()->origin().i() + 1;
  j = j - boxPtr()->origin().j() + 1;
  k = k - boxPtr()->origin().k();
  auto projection{_projection_x_half(i) + _projection_y_half(j) +
                  _projection_z_int(k) - 0.5};
  auto index{static_cast<std::size_t>(projection)};
  auto weight{projection - index};
  return (1 - weight) * _hz_inc(index) + weight * _hz_inc(index + 1);
}

double TFSF::cax() {
  return calculationParamPtr()->timeParam()->dt() /
         (constant::EPSILON_0 * gridSpacePtr()->basedDx());
}

double TFSF::cay() {
  return calculationParamPtr()->timeParam()->dt() /
         (constant::EPSILON_0 * gridSpacePtr()->basedDy());
}

double TFSF::caz() {
  return calculationParamPtr()->timeParam()->dt() /
         (constant::EPSILON_0 * gridSpacePtr()->basedDz());
}

double TFSF::cbx() {
  return calculationParamPtr()->timeParam()->dt() /
         (constant::MU_0 * gridSpacePtr()->basedDx());
}

double TFSF::cby() {
  return calculationParamPtr()->timeParam()->dt() /
         (constant::MU_0 * gridSpacePtr()->basedDy());
}

double TFSF::cbz() {
  return calculationParamPtr()->timeParam()->dt() /
         (constant::MU_0 * gridSpacePtr()->basedDz());
}

void TFSF::initTransform() {
  auto sin_theta{sinTheta()};
  auto cos_theta{cosTheta()};
  auto sin_phi{sinPhi()};
  auto cos_phi{cosPhi()};
  auto sin_psi{sinPsi()};
  auto cos_psi{cosPsi()};

  _k = Vector{sin_theta * cos_phi, sin_theta * sin_phi, cos_theta};

  _rotation_matrix = xt::xtensor<double, 2>{
      {-_sin_phi, _cos_theta * _cos_phi, _sin_theta * _cos_phi},
      {_cos_phi, _cos_theta * _sin_phi, _sin_theta * _sin_phi},
      {0, -_sin_theta, _cos_theta}};

  _k_e = Vector{sin_psi, cos_psi, 0};
  _transform_e = xt::linalg::dot(_rotation_matrix, _k_e.data());
  _transform_h = xt::linalg::cross(_k.data(), _transform_e);
}

void TFSF::calculateProjection() {
  auto get_int{[](std::size_t n, double inc_point, double tfsf_origin, double k,
                  double ratio_delta) {
    xt::xarray<double> arr = xt::zeros<double>({n});
    for (std::size_t i = 0; i < n; ++i) {
      arr(i) = (i + tfsf_origin - inc_point) * k * ratio_delta;
    }
    return arr;
  }};
  auto get_half{[](std::size_t n, double tfsf_origin, double inc_point,
                   double k, double ratio_delta) {
    xt::xarray<double> arr = xt::zeros<double>({n});
    for (std::size_t i = 0; i < n; ++i) {
      arr(i) = (i + tfsf_origin - inc_point - 0.5) * k * ratio_delta;
    }
    return arr;
  }};

  auto nx{boxPtr()->size().i()};
  auto ny{boxPtr()->size().j()};
  auto nz{boxPtr()->size().k()};

  auto tfsf_origin_x{boxPtr()->origin().i()};
  auto tfsf_origin_y{boxPtr()->origin().j()};
  auto tfsf_origin_z{boxPtr()->origin().k()};

  auto temp_point{calculateInjectPostion()};
  auto extra_dis{2 * _k.data() / _ratio_delta};
  auto inc_point_x{temp_point.i() - extra_dis(0)};
  auto inc_point_y{temp_point.j() - extra_dis(1)};
  auto inc_point_z{temp_point.k() - extra_dis(2)};

  _projection_x_int =
      get_int(nx + 1, inc_point_x, tfsf_origin_x, k().x(), _ratio_delta);
  _projection_y_int =
      get_int(ny + 1, inc_point_y, tfsf_origin_y, k().y(), _ratio_delta);
  _projection_z_int =
      get_int(nz + 1, inc_point_z, tfsf_origin_z, k().z(), _ratio_delta);
  _projection_x_half =
      get_half(nx + 2, tfsf_origin_x, inc_point_x, k().x(), _ratio_delta);
  _projection_y_half =
      get_half(ny + 2, tfsf_origin_y, inc_point_y, k().y(), _ratio_delta);
  _projection_z_half =
      get_half(nz + 2, tfsf_origin_z, inc_point_z, k().z(), _ratio_delta);
}

Grid TFSF::calculateInjectPostion() {
  if (0 <= _k.x() && 0 <= _k.y() && 0 <= _k.z()) {
    return _box->origin();
  }

  if (0 <= _k.x() && 0 <= _k.y() && _k.z() < 0) {
    return {_box->origin().i(), _box->origin().j(), _box->end().k()};
  }

  if (0 <= _k.x() && _k.y() < 0 && 0 <= _k.z()) {
    return {_box->origin().i(), _box->end().j(), _box->origin().k()};
  }

  if (0 <= _k.x() && _k.y() < 0 && _k.z() < 0) {
    return {_box->origin().i(), _box->end().j(), _box->end().k()};
  }

  if (_k.x() < 0 && 0 <= _k.y() && 0 <= _k.z()) {
    return {_box->end().i(), _box->origin().j(), _box->origin().k()};
  }

  if (_k.x() < 0 && 0 <= _k.y() && _k.z() < 0) {
    return {_box->end().i(), _box->origin().j(), _box->end().k()};
  }

  if (_k.x() < 0 && _k.y() < 0 && 0 <= _k.z()) {
    return {_box->end().i(), _box->end().j(), _box->origin().k()};
  }

  if (_k.x() < 0 && _k.y() < 0 && _k.z() < 0) {
    return _box->end();
  }

  throw std::runtime_error("TFSF::calculateInjectPostion() failed");
}

}  // namespace xfdtd
