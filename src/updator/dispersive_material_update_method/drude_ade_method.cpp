#include "updator/dispersive_material_update_method/drude_ade_method.h"

#include <xfdtd/common/constant.h>
#include <xfdtd/common/type_define.h>
#include <xfdtd/material/dispersive_material.h>

#include <utility>

#include "material/dispersive_material_equation.h"
#include "updator/dispersive_material_update_method/dispersive_material_update_method.h"
#include "updator/update_scheme.h"

namespace xfdtd {

auto DrudeADEMethod::k(const DrudeEqDecision& eq, Real dt) -> Array1D<Real> {
  const auto& gamma = eq.gamma();
  return (1 - gamma * dt / 2) / (1 + gamma * dt / 2);
}

auto DrudeADEMethod::beta(const DrudeEqDecision& eq, Real dt) -> Array1D<Real> {
  const auto& omega_p = eq.omegaP();
  const auto& gamma = eq.gamma();
  return (constant::EPSILON_0 * omega_p * omega_p * dt / 2) /
         (1 + gamma * dt / 2);
}

DrudeADEMethod::DrudeADEMethod(Real epsilon_inf,
                               std::shared_ptr<DrudeEqDecision> drude_eq)
    : LinearDispersiveMaterialUpdateMethod{epsilon_inf},
      _drude_eq{std::move(drude_eq)} {}

DrudeADEMethod::~DrudeADEMethod() = default;

auto DrudeADEMethod::init(Real dt) -> void {
  auto&& k = this->k(*_drude_eq, dt);
  auto&& beta = this->beta(*_drude_eq, dt);
  _sum_beta = std::accumulate(beta.begin(), beta.end(), 0.0);

  _coeff_j_j = k;
  _coeff_j_e = beta;
  _coeff_j_sum_j = 0.5 * (1 + _coeff_j_j);
}

auto DrudeADEMethod::correctCoeff(
    Index i, Index j, Index k, const GridSpace* grid_space,
    CalculationParam* calculation_param) const -> void {
  const auto& grid = grid_space->gridWithMaterial()(i, j, k);
  const auto sigma_e = calculation_param->materialParam()
                           ->materialArray()
                           .at(grid.materialIndex())
                           ->emProperty()
                           .sigmaE();
  const auto epsilon_inf = _epsilon_inf;
  const auto dt = calculation_param->timeParam()->dt();

  const auto a =
      this->a(epsilon_inf, constant::EPSILON_0, _sum_beta, dt, sigma_e);

  const auto b =
      this->b(epsilon_inf, constant::EPSILON_0, _sum_beta, dt, sigma_e);

  auto correct_func = [&a, &b](const auto da, const auto db, auto& cece,
                               auto& cecha, auto& cechb) {
    cece = a;
    cecha = -b / db;
    cechb = b / da;
  };

  const auto dx = grid_space->hSizeX()(i);
  const auto dy = grid_space->hSizeY()(j);
  const auto dz = grid_space->hSizeZ()(k);

  auto& cexe = calculation_param->fdtdCoefficient()->cexe()(i, j, k);
  auto& cexhy = calculation_param->fdtdCoefficient()->cexhy().at(i, j, k);
  auto& cexhz = calculation_param->fdtdCoefficient()->cexhz().at(i, j, k);

  correct_func(dy, dz, cexe, cexhy, cexhz);

  auto& ceye = calculation_param->fdtdCoefficient()->ceye()(i, j, k);
  auto& ceyhz = calculation_param->fdtdCoefficient()->ceyhz().at(i, j, k);
  auto& ceyhx = calculation_param->fdtdCoefficient()->ceyhx().at(i, j, k);

  correct_func(dz, dx, ceye, ceyhz, ceyhx);

  auto& ceze = calculation_param->fdtdCoefficient()->ceze()(i, j, k);
  auto& cezhx = calculation_param->fdtdCoefficient()->cezhx().at(i, j, k);
  auto& cezhy = calculation_param->fdtdCoefficient()->cezhy().at(i, j, k);

  correct_func(dx, dy, ceze, cezhx, cezhy);
}

auto DrudeADEMethod::initUpdate(
    const GridSpace* grid_space,
    std::shared_ptr<const CalculationParam> calculation_param,
    std::shared_ptr<EMF> emf, Index m_index, const IndexTask& task) -> void {
  if (!task.valid()) {
    throw XFDTDLinearDispersiveMaterialEquationException(
        "IndexTask is not valid.");
  }

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

  const auto b =
      this->b(_epsilon_inf, constant::EPSILON_0, _sum_beta, dt, sigma_e);

  _coeff_e_j_sum = -b;

  const auto num_p = _drude_eq->numPoles();

  const auto nx = task.xRange().size();
  const auto ny = task.yRange().size();
  const auto nz = task.zRange().size();

  _jx = xt::zeros<Real>({nx, ny + 1, nz + 1, num_p});
  _jy = xt::zeros<Real>({nx + 1, ny, nz + 1, num_p});
  _jz = xt::zeros<Real>({nx + 1, ny + 1, nz, num_p});
}

auto DrudeADEMethod::updateEx(Index i, Index j, Index k) -> void {
  auto& ex{_emf->ex()};
  const auto& cexe{_calculation_param->fdtdCoefficient()->cexe()};
  const auto& cexhy{_calculation_param->fdtdCoefficient()->cexhy()};
  const auto& cexhz{_calculation_param->fdtdCoefficient()->cexhz()};
  const auto& hy{_emf->hy()};
  const auto& hz{_emf->hz()};

  auto j_sum = calculateJSum(i, j, k, _jx);
  const auto e_cur = ex(i, j, k);
  auto&& e_next = ex(i, j, k);

  e_next =
      eNext(cexe(i, j, k), ex(i, j, k), cexhy(i, j, k), hy(i, j, k),
            hy(i, j, k - 1), cexhz(i, j, k), hz(i, j, k), hz(i, j - 1, k)) +
      _coeff_e_j_sum * j_sum;

  updateJ(i, j, k, e_next, e_cur, _jx);
}

auto DrudeADEMethod::updateEy(Index i, Index j, Index k) -> void {
  auto& ey{_emf->ey()};
  const auto& ceye{_calculation_param->fdtdCoefficient()->ceye()};
  const auto& ceyhz{_calculation_param->fdtdCoefficient()->ceyhz()};
  const auto& ceyhx{_calculation_param->fdtdCoefficient()->ceyhx()};
  const auto& hz{_emf->hz()};
  const auto& hx{_emf->hx()};

  auto j_sum = calculateJSum(i, j, k, _jy);
  const auto e_cur = ey(i, j, k);
  auto&& e_next = ey(i, j, k);

  e_next =
      eNext(ceye(i, j, k), ey(i, j, k), ceyhz(i, j, k), hz(i, j, k),
            hz(i - 1, j, k), ceyhx(i, j, k), hx(i, j, k), hx(i, j, k - 1)) +
      _coeff_e_j_sum * j_sum;

  updateJ(i, j, k, e_next, e_cur, _jy);
}

auto DrudeADEMethod::updateEz(Index i, Index j, Index k) -> void {
  auto& ez{_emf->ez()};
  const auto& ceze{_calculation_param->fdtdCoefficient()->ceze()};
  const auto& cezhx{_calculation_param->fdtdCoefficient()->cezhx()};
  const auto& cezhy{_calculation_param->fdtdCoefficient()->cezhy()};
  const auto& hx{_emf->hx()};
  const auto& hy{_emf->hy()};

  auto j_sum = calculateJSum(i, j, k, _jz);
  const auto e_cur = ez(i, j, k);
  auto&& e_next = ez(i, j, k);

  e_next =
      eNext(ceze(i, j, k), ez(i, j, k), cezhx(i, j, k), hx(i, j, k),
            hx(i, j - 1, k), cezhy(i, j, k), hy(i, j, k), hy(i - 1, j, k)) +
      _coeff_e_j_sum * j_sum;

  updateJ(i, j, k, e_next, e_cur, _jz);
}

auto DrudeADEMethod::updateTEM(Index i, Index j, Index k) -> void {
  auto& ex{_emf->ex()};
  const auto& cexe{_calculation_param->fdtdCoefficient()->cexe()};
  const auto& cexhy{_calculation_param->fdtdCoefficient()->cexhy()};
  const auto& hy{_emf->hy()};

  auto j_sum = calculateJSum(i, j, k, _jx);
  const auto e_cur = ex(i, j, k);
  auto&& e_next = ex(i, j, k);

  e_next = eNext(cexe(i, j, k), ex(i, j, k), cexhy(i, j, k), hy(i, j, k),
                 hy(i, j, k - 1), 0.0, 0.0, 0.0) +
           _coeff_e_j_sum * j_sum;

  updateJ(i, j, k, e_next, e_cur, _jx);
}

auto DrudeADEMethod::calculateJSum(Index i, Index j, Index k,
                                   const Array4D<Real>& j_arr) const -> Real {
  Real j_sum{0};
  auto num_p = j_arr.shape().back();

  auto local_i = i - _is;
  auto local_j = j - _js;
  auto local_k = k - _ks;

  for (Index p{0}; p < num_p; ++p) {
    j_sum += _coeff_j_sum_j(p) * j_arr(local_i, local_j, local_k, p);
  }

  return j_sum;
}

auto DrudeADEMethod::updateJ(Index i, Index j, Index k, Real e_next, Real e_cur,
                             Array4D<Real>& j_arr) -> void {
  auto num_p = j_arr.shape().back();
  auto local_i = i - _is;
  auto local_j = j - _js;
  auto local_k = k - _ks;

  for (Index p{0}; p < num_p; ++p) {
    j_arr(local_i, local_j, local_k, p) =
        _coeff_j_j(p) * j_arr(local_i, local_j, local_k, p) +
        _coeff_j_e(p) * (e_next + e_cur);
  }
}

}  // namespace xfdtd
