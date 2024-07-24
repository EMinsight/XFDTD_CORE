#include <xfdtd/common/constant.h>
#include <xfdtd/common/index_task.h>
#include <xfdtd/grid_space/grid_space.h>
#include <xfdtd/waveform_source/tfsf.h>

#include <cmath>
#include <cstdlib>

#include "xfdtd/common/type_define.h"
#include "xfdtd/coordinate_system/coordinate_system.h"

namespace xfdtd {

TFSF::TFSF(std::size_t x, std::size_t y, std::size_t z, Real theta, Real phi,
           Real psi, std::unique_ptr<Waveform> waveform)
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

void TFSF::initTimeDependentVariable() {
  waveform()->init(calculationParamPtr()->timeParam()->eTime() -
                   calculationParamPtr()->timeParam()->dt());

  auto fdtd_1d = [](auto&& e, auto&& h, const auto ceie, const auto ceih,
                    const auto chih, const auto chie, const auto abc_coeff_0,
                    const auto abc_coeff_1, const auto& source) {
    const auto nt = h.shape()[0];
    const auto nl = h.shape()[1];

    e(0, 0) = source(0);
    for (Index l = 1; l < nl; ++l) {
      h(0, l) = chih * h(0, l) + chie * (e(0, l + 1) - e(0, l));
    }

    for (Index t = 1; t < nt; ++t) {
      for (Index l = 1; l < nl; ++l) {
        e(t, 0) = source(t);
        e(t, l) = ceie * e(t - 1, l) + ceih * (h(t - 1, l) - h(t - 1, l - 1));
      }

      // abc
      Real e_l_prev_t_prev_prev = (t < 2) ? 0.0 : e(t - 2, nl - 1);
      Real e_l_t_prev_prev = (t < 2) ? 0.0 : e(t - 2, nl);
      e(t, nl) = -e_l_prev_t_prev_prev +
                 abc_coeff_0 * (e(t, nl - 1) + e_l_t_prev_prev) +
                 abc_coeff_1 * (e(t - 1, nl) + e(t - 1, nl - 1));

      for (Index l = 0; l < nl; ++l) {
        h(t, l) = chih * h(t - 1, l) + chie * (e(t, l + 1) - e(t, l));
      }
    }
  };

  _e_inc = xt::zeros<Real>(
      {calculationParamPtr()->timeParam()->size(), _auxiliary_size});

  _h_inc = xt::zeros<Real>(
      {calculationParamPtr()->timeParam()->size(), _auxiliary_size - 1});

  fdtd_1d(_e_inc, _h_inc, _ceie, _ceihi, _chih, _chiei, _abc_coff_0,
          _abc_coff_1, waveform()->value());
}

std::size_t TFSF::x() const { return _x; }

std::size_t TFSF::y() const { return _y; }

std::size_t TFSF::z() const { return _z; }

Real TFSF::theta() const { return _theta; }

Real TFSF::phi() const { return _phi; }

Real TFSF::psi() const { return _psi; }

Real TFSF::sinTheta() const { return _sin_theta; }

Real TFSF::cosTheta() const { return _cos_theta; }

Real TFSF::sinPhi() const { return _sin_phi; }

Real TFSF::cosPhi() const { return _cos_phi; }

Real TFSF::sinPsi() const { return _sin_psi; }

Real TFSF::cosPsi() const { return _cos_psi; }

Vector TFSF::k() const { return _k; }

GridBox TFSF::globalBox() const { return _global_box; }

auto TFSF::globalTask() const -> IndexTask {
  const auto is = globalBox().origin().i();
  const auto ie = globalBox().end().i();
  const auto js = globalBox().origin().j();
  const auto je = globalBox().end().j();
  const auto ks = globalBox().origin().k();
  const auto ke = globalBox().end().k();
  return makeIndexTask(makeRange(is, ie), makeRange(js, je), makeRange(ks, ke));
}

void TFSF::defaultInit(std::shared_ptr<GridSpace> grid_space,
                       std::shared_ptr<CalculationParam> calculation_param,
                       std::shared_ptr<EMF> emf) {
  WaveformSource::defaultInit(std::move(grid_space),
                              std::move(calculation_param), std::move(emf));
  if (gridSpacePtr()->type() != GridSpace::Type::UNIFORM) {
    throw std::runtime_error("TFSF only supports uniform grid space");
  }

  auto global_grid_space = gridSpacePtr()->globalGridSpace();
  if (global_grid_space == nullptr) {
    throw std::runtime_error("TFSF need global grid space");
  }

  auto origin_x{global_grid_space->box().origin().i() + x()};
  auto origin_y{global_grid_space->box().origin().j() + y()};
  auto origin_z{global_grid_space->box().origin().k() + z()};
  auto size_x{global_grid_space->box().size().i() - 2 * x()};
  auto size_y{global_grid_space->box().size().j() - 2 * y()};
  auto size_z{global_grid_space->box().size().k() - 2 * z()};

  _global_box =
      GridBox{Grid{origin_x, origin_y, origin_z}, Grid{size_x, size_y, size_z}};

  _ratio_delta =
      1 / std::sqrt(std::pow(sinTheta(), 4) *
                        (std::pow(cosPhi(), 4) + std::pow(sinPhi(), 4)) +
                    std::pow(cosTheta(), 4));
  _auxiliary_size =
      std::ceil(_ratio_delta *
                (std::sqrt(pow(size_x, 2) + pow(size_y, 2) + pow(size_z, 2)))) +
      4 + 1;

  calculateProjection();

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

Real TFSF::cax() const {
  return calculationParam()->timeParam()->dt() /
         (constant::EPSILON_0 * gridSpace()->basedDx());
}

Real TFSF::cay() const {
  return calculationParam()->timeParam()->dt() /
         (constant::EPSILON_0 * gridSpace()->basedDy());
}

Real TFSF::caz() const {
  return calculationParam()->timeParam()->dt() /
         (constant::EPSILON_0 * gridSpace()->basedDz());
}

Real TFSF::cbx() const {
  return calculationParam()->timeParam()->dt() /
         (constant::MU_0 * gridSpace()->basedDx());
}

Real TFSF::cby() const {
  return calculationParam()->timeParam()->dt() /
         (constant::MU_0 * gridSpace()->basedDy());
}

Real TFSF::cbz() const {
  return calculationParam()->timeParam()->dt() /
         (constant::MU_0 * gridSpace()->basedDz());
}

auto TFSF::nodeTask() const -> IndexTask {
  auto node_global_task = nodeGlobalTask();
  if (!node_global_task.valid()) {
    return {};
  }

  auto g_box = gridSpace()->globalBox();
  auto offset = g_box.origin();
  auto node_task = makeIndexTask(node_global_task.xRange() - offset.i(),
                                 node_global_task.yRange() - offset.j(),
                                 node_global_task.zRange() - offset.k());
  return node_task;
}

auto TFSF::nodeGlobalTask() const -> IndexTask {
  auto g_box = gridSpace()->globalBox();
  auto g_task = makeTask(makeIndexRange(g_box.origin().i(), g_box.end().i()),
                         makeIndexRange(g_box.origin().j(), g_box.end().j()),
                         makeIndexRange(g_box.origin().k(), g_box.end().k()));
  auto tfsf_g_task = globalTask();
  auto intersection = taskIntersection(g_task, tfsf_g_task);
  return intersection.has_value() ? intersection.value() : IndexTask{};
}

void TFSF::initTransform() {
  auto sin_theta{sinTheta()};
  auto cos_theta{cosTheta()};
  auto sin_phi{sinPhi()};
  auto cos_phi{cosPhi()};
  auto sin_psi{sinPsi()};
  auto cos_psi{cosPsi()};

  _k = Vector{sin_theta * cos_phi, sin_theta * sin_phi, cos_theta};

  auto a = Vector{-_sin_phi, _cos_theta * _cos_phi, _sin_theta * _cos_phi};
  auto b = Vector{_cos_phi, _cos_theta * _sin_phi, _sin_theta * _sin_phi};
  auto c = Vector{0, -_sin_theta, _cos_theta};

  _k_e = Vector{sin_psi, cos_psi, 0};

  _transform_e = Vector{
      a.dot(_k_e),
      b.dot(_k_e),
      c.dot(_k_e),
  };

  _transform_h = _k.cross(_transform_e);
}

void TFSF::calculateProjection() {
  auto get_int{[](std::size_t n, Real inc_point, Real tfsf_origin, Real k,
                  Real ratio_delta) {
    Array1D<Real> arr = xt::zeros<Real>({n});
    for (std::size_t i = 0; i < n; ++i) {
      arr(i) = (i + tfsf_origin - inc_point) * k * ratio_delta;
    }
    return arr;
  }};
  auto get_half{[](std::size_t n, Real tfsf_origin, Real inc_point, Real k,
                   Real ratio_delta) {
    Array1D<Real> arr = xt::zeros<Real>({n});
    for (std::size_t i = 0; i < n; ++i) {
      arr(i) = (i + tfsf_origin - inc_point - 0.5) * k * ratio_delta;
    }
    return arr;
  }};

  const auto& box = globalBox();

  auto nx{box.size().i()};
  auto ny{box.size().j()};
  auto nz{box.size().k()};

  auto tfsf_origin_x{box.origin().i()};
  auto tfsf_origin_y{box.origin().j()};
  auto tfsf_origin_z{box.origin().k()};

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
  const auto& box = globalBox();
  if (0 <= _k.x() && 0 <= _k.y() && 0 <= _k.z()) {
    return box.origin();
  }

  if (0 <= _k.x() && 0 <= _k.y() && _k.z() < 0) {
    return {box.origin().i(), box.origin().j(), box.end().k()};
  }

  if (0 <= _k.x() && _k.y() < 0 && 0 <= _k.z()) {
    return {box.origin().i(), box.end().j(), box.origin().k()};
  }

  if (0 <= _k.x() && _k.y() < 0 && _k.z() < 0) {
    return {box.origin().i(), box.end().j(), box.end().k()};
  }

  if (_k.x() < 0 && 0 <= _k.y() && 0 <= _k.z()) {
    return {box.end().i(), box.origin().j(), box.origin().k()};
  }

  if (_k.x() < 0 && 0 <= _k.y() && _k.z() < 0) {
    return {box.end().i(), box.origin().j(), box.end().k()};
  }

  if (_k.x() < 0 && _k.y() < 0 && 0 <= _k.z()) {
    return {box.end().i(), box.end().j(), box.origin().k()};
  }

  if (_k.x() < 0 && _k.y() < 0 && _k.z() < 0) {
    return box.end();
  }

  throw std::runtime_error("TFSF::calculateInjectPostion() failed");
}

}  // namespace xfdtd
