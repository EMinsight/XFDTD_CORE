#include "updator/dispersive_material_update_method/m_lorentz_ade_method.h"

#include <xtensor.hpp>

#include "material/dispersive_material_equation.h"
#include "updator/update_scheme.h"
#include "updator/updator.h"

namespace xfdtd {

auto MLorentzADEMethod::a(const MLorentzEqDecision& m_lorentz_equation, Real dt)
    -> Array1D<Real> {
  auto&& a_1 = m_lorentz_equation.a1() * constant::EPSILON_0;
  auto&& a_0 = m_lorentz_equation.a0() * constant::EPSILON_0;
  return ((a_1 / (dt * dt)) + (0.5 * a_0 / dt));
}

auto MLorentzADEMethod::b(const MLorentzEqDecision& m_lorentz_equation, Real dt)
    -> Array1D<Real> {
  auto&& a_1 = m_lorentz_equation.a1() * constant::EPSILON_0;
  return (-2 * a_1) / (dt * dt);
}

auto MLorentzADEMethod::c(const MLorentzEqDecision& m_lorentz_equation, Real dt)
    -> Array1D<Real> {
  auto&& a_1 = m_lorentz_equation.a1() * constant::EPSILON_0;
  auto&& a_0 = m_lorentz_equation.a0() * constant::EPSILON_0;
  return ((a_1 / (dt * dt)) - (0.5 * a_0 / dt));
}

auto MLorentzADEMethod::d(const MLorentzEqDecision& m_lorentz_equation, Real dt)
    -> Array1D<Real> {
  auto&& b_2 = m_lorentz_equation.b2();
  auto&& b_0 = m_lorentz_equation.b0();
  return ((2 * b_2 / (dt * dt)) - b_0);
}

auto MLorentzADEMethod::e(const MLorentzEqDecision& m_lorentz_equation, Real dt)
    -> Array1D<Real> {
  auto&& b_2 = m_lorentz_equation.b2();
  auto&& b_1 = m_lorentz_equation.b1();
  return ((-b_2 / (dt * dt)) + (b_1 / (2 * dt)));
}

auto MLorentzADEMethod::f(const MLorentzEqDecision& m_lorentz_equation, Real dt)
    -> Array1D<Real> {
  auto&& b_2 = m_lorentz_equation.b2();
  auto&& b_1 = m_lorentz_equation.b1();
  return ((b_2 / (dt * dt)) + (b_1 / (2 * dt)));
}

MLorentzADEMethod::MLorentzADEMethod(
    Real epsilon_inf, std::shared_ptr<MLorentzEqDecision> m_lorentz_equation)
    : LinearDispersiveMaterialUpdateMethod{epsilon_inf},
      _m_lorentz_equation{std::move(m_lorentz_equation)} {
  if (_m_lorentz_equation == nullptr) {
    throw XFDTDUpdatorException("MLorentzEqDecision is nullptr.");
  }
}

auto MLorentzADEMethod::init(Real dt) -> void {
  _a = a(*_m_lorentz_equation, dt);
  _b = b(*_m_lorentz_equation, dt);
  _c = c(*_m_lorentz_equation, dt);
  _d = d(*_m_lorentz_equation, dt);
  _e = e(*_m_lorentz_equation, dt);
  _f = f(*_m_lorentz_equation, dt);
  _coeff_j_e_n = _a / _f;
  _coeff_j_e_c = _b / _f;
  _coeff_j_e_p = _c / _f;
  _coeff_j_j_c = _d / _f;
  _coeff_j_j_p = _e / _f;
  _coeff_j_sum_j_c = 0.5 * (1 + _d / _f);
  _coeff_j_sum_j_p = 0.5 * _e / _f;
  auto omega_p = xt::sqrt(_m_lorentz_equation->b0());
  auto delta_epsilon = _m_lorentz_equation->a0() / _m_lorentz_equation->b0();
  auto nv = 0.5 * _m_lorentz_equation->b1();
}

auto MLorentzADEMethod::correctCoeff(Index i, Index j, Index k,
                                     const GridSpace* grid_space,
                                     CalculationParam* calculation_param) const
    -> void {
  Real sum_a_f = 0;
  for (Index p{0}; p < _a.size(); ++p) {
    sum_a_f += _a(p) / _f(p);
  }
  Real sum_b_f = 0;
  for (Index p{0}; p < _b.size(); ++p) {
    sum_b_f += _b(p) / _f(p);
  }

  auto correct_func = [epsilon_inf = _epsilon_inf, sum_a_f, sum_b_f](
                          const auto& dt, const auto& da, const auto& db,
                          const auto& sigma_e, auto& cece, auto& cecha,
                          auto& cechb) {
    const auto epsilon_0 = xfdtd::constant::EPSILON_0;
    auto&& temp =
        (epsilon_inf * epsilon_0) / dt + 0.5 * sum_a_f + 0.5 * sigma_e;
    cece =
        ((epsilon_inf * epsilon_0) / dt - 0.5 * sum_b_f - 0.5 * sigma_e) / temp;
    cecha = -1 / (temp * db);
    cechb = 1 / (temp * da);
  };

  const auto dt = calculation_param->timeParam()->dt();
  const auto& grid = grid_space->gridWithMaterial()(i, j, k);
  const auto sigma_e = calculation_param->materialParam()
                           ->materialArray()
                           .at(grid->materialIndex())
                           ->emProperty()
                           .sigmaE();
  const auto dx = grid_space->hSizeX().at(i);
  const auto dy = grid_space->hSizeY().at(j);
  const auto dz = grid_space->hSizeZ().at(k);

  auto& cexe = calculation_param->fdtdCoefficient()->cexe().at(i, j, k);
  auto& cexhy = calculation_param->fdtdCoefficient()->cexhy().at(i, j, k);
  auto& cexhz = calculation_param->fdtdCoefficient()->cexhz().at(i, j, k);

  correct_func(dt, dy, dz, sigma_e, cexe, cexhy, cexhz);

  auto& ceye = calculation_param->fdtdCoefficient()->ceye().at(i, j, k);
  auto& ceyhz = calculation_param->fdtdCoefficient()->ceyhz().at(i, j, k);
  auto& ceyhx = calculation_param->fdtdCoefficient()->ceyhx().at(i, j, k);

  correct_func(dt, dz, dx, sigma_e, ceye, ceyhz, ceyhx);

  auto& ceze = calculation_param->fdtdCoefficient()->ceze().at(i, j, k);
  auto& cezhx = calculation_param->fdtdCoefficient()->cezhx().at(i, j, k);
  auto& cezhy = calculation_param->fdtdCoefficient()->cezhy().at(i, j, k);

  correct_func(dt, dx, dy, sigma_e, ceze, cezhx, cezhy);
}

auto MLorentzADEMethod::initUpdate(
    const GridSpace* grid_space,
    std::shared_ptr<const CalculationParam> calculation_param,
    std::shared_ptr<EMF> emf, Index m_index, const IndexTask& task) -> void {
  if (!task.valid()) {
    throw XFDTDUpdatorException("IndexTask is not valid.");
  }

  checkStabilityCondition(calculation_param->timeParam()->dt());

  _is = task.xRange().start();
  _js = task.yRange().start();
  _ks = task.zRange().start();

  _calculation_param = std::move(calculation_param);
  _emf = std::move(emf);

  const auto dt = _calculation_param->timeParam()->dt();

  const auto sigma_e = _calculation_param->materialParam()
                           ->materialArray()
                           .at(m_index)
                           ->emProperty()
                           .sigmaE();
  Real sum_a_f = 0;
  for (Index p{0}; p < _f.size(); ++p) {
    sum_a_f += _a(p) / _f(p);
  }
  Real sum_c_f = 0;
  for (Index p{0}; p < _f.size(); ++p) {
    sum_c_f += _c(p) / _f(p);
  }
  Real temp =
      (_epsilon_inf * constant::EPSILON_0) / dt + 0.5 * sum_a_f + 0.5 * sigma_e;

  _coeff_e_j_sum = -1 / temp;
  _coeff_e_e_p = -0.5 * sum_c_f / temp;

  const Material* m_ptr =
      _calculation_param->materialParam()->materialArray().at(m_index).get();

  auto linear_dispersive_material_prt =
      dynamic_cast<const LinearDispersiveMaterial*>(m_ptr);
  if (linear_dispersive_material_prt == nullptr) {
    throw XFDTDUpdatorException("Material is not LinearDispersiveMaterial.");
  }

  auto num_p = linear_dispersive_material_prt->numPoles();

  const auto nx = task.xRange().size();
  const auto ny = task.yRange().size();
  const auto nz = task.zRange().size();

  _ex_prev = xt::zeros<Real>({nx, ny + 1, nz + 1});
  _ey_prev = xt::zeros<Real>({nx + 1, ny, nz + 1});
  _ez_prev = xt::zeros<Real>({nx + 1, ny + 1, nz});

  _jx = xt::zeros<Real>({nx, ny + 1, nz + 1, num_p});
  _jy = xt::zeros<Real>({nx + 1, ny, nz + 1, num_p});
  _jz = xt::zeros<Real>({nx + 1, ny + 1, nz, num_p});
  _jx_prev = xt::zeros_like(_jx);
  _jy_prev = xt::zeros_like(_jy);
  _jz_prev = xt::zeros_like(_jz);
}

auto MLorentzADEMethod::updateEx(Index i, Index j, Index k) -> void {
  auto& ex = _emf->ex();
  const auto& cexe = _calculation_param->fdtdCoefficient()->cexe();
  const auto& cexhy = _calculation_param->fdtdCoefficient()->cexhy();
  const auto& cexhz = _calculation_param->fdtdCoefficient()->cexhz();
  const auto& hy = _emf->hy();
  const auto& hz = _emf->hz();

  auto j_sum = calculateJSum(i, j, k, _jx, _jx_prev);

  auto local_i = i - _is;
  auto local_j = j - _js;
  auto local_k = k - _ks;
  const auto e_cur = ex(i, j, k);
  auto& e_prev = _ex_prev(local_i, local_j, local_k);

  auto& e_next = ex(i, j, k);

  e_next =
      _coeff_e_e_p * e_prev +
      eNext(cexe(i, j, k), ex(i, j, k), cexhy(i, j, k), hy(i, j, k),
            hy(i, j, k - 1), cexhz(i, j, k), hz(i, j, k), hz(i, j - 1, k)) +
      _coeff_e_j_sum * j_sum;

  updateJ(i, j, k, e_next, e_cur, e_prev, _jx, _jx_prev);

  e_prev = e_cur;
}

auto MLorentzADEMethod::updateEy(Index i, Index j, Index k) -> void {
  auto& ey = _emf->ey();
  const auto& ceye = _calculation_param->fdtdCoefficient()->ceye();
  const auto& ceyhz = _calculation_param->fdtdCoefficient()->ceyhz();
  const auto& ceyhx = _calculation_param->fdtdCoefficient()->ceyhx();
  const auto& hz = _emf->hz();
  const auto& hx = _emf->hx();

  auto j_sum = calculateJSum(i, j, k, _jy, _jy_prev);

  auto local_i = i - _is;
  auto local_j = j - _js;
  auto local_k = k - _ks;
  const auto e_cur = ey(i, j, k);
  auto& e_prev = _ey_prev(local_i, local_j, local_k);

  auto& e_next = ey(i, j, k);

  e_next =
      _coeff_e_e_p * e_prev +
      eNext(ceye(i, j, k), ey(i, j, k), ceyhz(i, j, k), hz(i, j, k),
            hz(i - 1, j, k), ceyhx(i, j, k), hx(i, j, k), hx(i, j, k - 1)) +
      _coeff_e_j_sum * j_sum;

  updateJ(i, j, k, e_next, e_cur, e_prev, _jy, _jy_prev);

  e_prev = e_cur;
}

auto MLorentzADEMethod::updateEz(Index i, Index j, Index k) -> void {
  auto& ez = _emf->ez();
  const auto& ceze = _calculation_param->fdtdCoefficient()->ceze();
  const auto& cezhx = _calculation_param->fdtdCoefficient()->cezhx();
  const auto& cezhy = _calculation_param->fdtdCoefficient()->cezhy();
  const auto& hx = _emf->hx();
  const auto& hy = _emf->hy();

  auto j_sum = calculateJSum(i, j, k, _jz, _jz_prev);

  auto local_i = i - _is;
  auto local_j = j - _js;
  auto local_k = k - _ks;
  const auto e_cur = ez(i, j, k);
  auto& e_prev = _ez_prev(local_i, local_j, local_k);

  auto& e_next = ez(i, j, k);

  e_next =
      _coeff_e_e_p * e_prev +
      eNext(ceze(i, j, k), ez(i, j, k), cezhx(i, j, k), hx(i, j, k),
            hx(i, j - 1, k), cezhy(i, j, k), hy(i, j, k), hy(i - 1, j, k)) +
      _coeff_e_j_sum * j_sum;

  updateJ(i, j, k, e_next, e_cur, e_prev, _jz, _jz_prev);

  e_prev = e_cur;
}

auto MLorentzADEMethod::updateTEM(Index i, Index j, Index k) -> void {
  auto& ex = _emf->ex();
  const auto& cexe = _calculation_param->fdtdCoefficient()->cexe();
  const auto& cexhy = _calculation_param->fdtdCoefficient()->cexhy();
  const auto& hy = _emf->hy();

  auto j_sum = calculateJSum(i, j, k, _jx, _jx_prev);

  auto local_i = i - _is;
  auto local_j = j - _js;
  auto local_k = k - _ks;
  const auto e_cur = ex(i, j, k);
  auto& e_prev = _ex_prev(local_i, local_j, local_k);

  auto& e_next = ex(i, j, k);

  e_next = _coeff_e_e_p * e_prev +
           eNext(cexe(i, j, k), ex(i, j, k), cexhy(i, j, k), hy(i, j, k),
                 hy(i, j, k - 1), 0.0, 0.0, 0.0) +
           _coeff_e_j_sum * j_sum;

  updateJ(i, j, k, e_next, e_cur, e_prev, _jx, _jx_prev);

  e_prev = e_cur;
}

auto MLorentzADEMethod::calculateJSum(Index i, Index j, Index k,
                                      const Array4D<Real>& j_arr,
                                      const Array4D<Real>& j_prev_arr) -> Real {
  auto num_p = j_arr.shape().back();
  auto local_i = i - _is;
  auto local_j = j - _js;
  auto local_k = k - _ks;

  Real sum = 0;
  for (Index p{0}; p < num_p; ++p) {
    sum += _coeff_j_sum_j_c(p) * j_arr(local_i, local_j, local_k, p) +
           _coeff_j_sum_j_p(p) * j_prev_arr(local_i, local_j, local_k, p);
  }

  return sum;
}

auto MLorentzADEMethod::updateJ(Index i, Index j, Index k, Real e_next,
                                Real e_cur, Real e_prev, Array4D<Real>& j_arr,
                                Array4D<Real>& j_prev_arr) -> void {
  auto num_p = j_arr.shape().back();
  auto local_i = i - _is;
  auto local_j = j - _js;
  auto local_k = k - _ks;

  for (Index p{0}; p < num_p; ++p) {
    auto& j_prev = j_prev_arr(local_i, local_j, local_k, p);
    const auto j_cur = j_arr(local_i, local_j, local_k, p);

    j_arr(local_i, local_j, local_k, p) =
        _coeff_j_e_n(p) * e_next + _coeff_j_e_c(p) * e_cur +
        _coeff_j_e_p(p) * e_prev + _coeff_j_j_c(p) * j_cur +
        _coeff_j_j_p(p) * j_prev;

    j_prev = j_cur;
  }
}

auto MLorentzADEMethod::checkStabilityCondition(Real dt) -> void {
  bool c_1 = !xt::any(_m_lorentz_equation->b0() < 0);
  bool c_2 = !xt::any(_m_lorentz_equation->b1() < 0);
  bool c_3 = !xt::any(_m_lorentz_equation->b2() < 0);
  // assume that v^2 is always satisfied.
  bool c_4 =
      !xt::any((_m_lorentz_equation->a0() * _m_lorentz_equation->b1() -
                _m_lorentz_equation->a1() * _m_lorentz_equation->b0()) < 0);
  bool c_5 = !xt::any((4 * _m_lorentz_equation->b2() -
                       _m_lorentz_equation->b0() * dt * dt) < 0);

  constexpr std::string_view err_c_1 =
      "b_0 must be greater than or equal to 0.";
  constexpr std::string_view err_c_2 =
      "b_1 must be greater than or equal to 0.";
  constexpr std::string_view err_c_3 =
      "b_2 must be greater than or equal to 0.";
  constexpr std::string_view err_c_4 = "Q must be greater than or equal to 0.";
  constexpr std::string_view err_c_5 =
      "Q_3 must be greater than or equal to 0.";

  std::stringstream ss;
  ss << "Check stability condition.\n";
  if (!c_1) {
    ss << "Warning: " << err_c_1 << "\n";
  }
  if (!c_2) {
    ss << "Warning: " << err_c_2 << "\n";
  }
  if (!c_3) {
    ss << "Warning: " << err_c_3 << "\n";
  }
  if (!c_4) {
    ss << "Warning: " << err_c_4 << "\n";
  }
  if (!c_5) {
    ss << "Warning: " << err_c_5 << "\n";
  }
  if (c_1 && c_2 && c_3 && c_4 && c_5) {
    ss << "All conditions are satisfied.\n";
  }

  std::cout << ss.str();
}

}  // namespace xfdtd