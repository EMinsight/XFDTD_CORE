#include <xfdtd/common/type_define.h>
#include <xfdtd/material/ade_method/drude_ade_method.h>
#include <xfdtd/material/dispersive_material_equation/dispersive_material_equation.h>

#include <sstream>
#include <xtensor/xbuilder.hpp>

namespace xfdtd {

class DrudeADEMethodCoefficient {
 public:
  static auto k(const Array1D<Real>& gamma, Real dt) -> Array1D<Real> {
    return (1 - gamma * dt / 2) / (1 + gamma * dt / 2);
  }

  static auto beta(const Array1D<Real>& omega_p, const Array1D<Real>& gamma,
                   Real dt) -> Array1D<Real> {
    return (constant::EPSILON_0 * omega_p * omega_p * dt / 2) /
           (1 + gamma * dt / 2);
  }

  static auto sumBeta(const Array1D<Real>& beta) -> Real {
    return std::accumulate(beta.begin(), beta.end(), Real{0.0});
  }

  static auto a(const Real& epsilon_inf, const Real& epsilon_0,
                const Real& sum_beta, const Real& dt,
                const Real& sigma) -> Real {
    return (2 * epsilon_0 * epsilon_inf - dt * sum_beta - dt * sigma) /
           (2 * epsilon_0 * epsilon_inf - dt * sum_beta + dt * sigma);
  }

  static auto b(const Real& epsilon_inf, const Real& epsilon_0,
                const Real& sum_beta, const Real& dt,
                const Real& sigma) -> Real {
    return (2 * dt) /
           (2 * epsilon_0 * epsilon_inf - dt * sum_beta + dt * sigma);
  }

  static auto cece(Real a) { return a; }

  static auto cecha(Real b, Real db) { return -b / db; }

  static auto chchb(Real b, Real da) { return b / da; }

  static auto coeffEJSum(Real b) { return -b; }
};

DrudeADEMethodStorage::DrudeADEMethodStorage(Index num_pole, Index nx, Index ny,
                                             Index nz)
    : ADEMethodStorage{num_pole} {
  _coeff_j_j = xt::zeros<Real>({nx, ny, nz, num_pole});
  _coeff_j_e = xt::zeros<Real>({nx, ny, nz, num_pole});
  _coeff_j_sum_j = xt::zeros<Real>({nx, ny, nz, num_pole});
  _coeff_e_j_sum = xt::zeros<Real>({nx, ny, nz});
  _jx_arr = xt::zeros<Real>({nx, ny, nz, num_pole});
  _jy_arr = xt::zeros<Real>({nx, ny, nz, num_pole});
  _jz_arr = xt::zeros<Real>({nx, ny, nz, num_pole});
}

auto DrudeADEMethodStorage::correctCoeff(
    Index i, Index j, Index k,
    const LinearDispersiveMaterial& linear_dispersive_material,
    const std::shared_ptr<const GridSpace>& grid_space,
    const std::shared_ptr<CalculationParam>& calculation_param) -> void {
  if (_num_pole < linear_dispersive_material.numPoles()) {
    std::stringstream ss;
    ss << "The number of poles is not enough. The number of poles is "
       << _num_pole << " but the number of poles in the material is "
       << linear_dispersive_material.numPoles();
    throw XFDTDLinearDispersiveMaterialException{ss.str()};
  }

  auto drude_eq = std::dynamic_pointer_cast<DrudeEqDecision>(
      linear_dispersive_material.equation());
  if (!drude_eq) {
    throw XFDTDLinearDispersiveMaterialException{
        "The equation is not Drude equation."};
  }

  const auto& omega_p = drude_eq->omegaP();
  const auto& gamma = drude_eq->gamma();
  const auto& epsilon_inf = linear_dispersive_material.epsilonInf();
  const auto& sigma = linear_dispersive_material.emProperty().sigmaE();

  const auto dt = calculation_param->timeParam()->dt();

  auto coeff_k = DrudeADEMethodCoefficient::k(gamma, dt);
  auto coeff_beta = DrudeADEMethodCoefficient::beta(omega_p, gamma, dt);
  const auto& sum_beta = DrudeADEMethodCoefficient::sumBeta(coeff_beta);

  for (Index p{0}; p < _num_pole; ++p) {
    _coeff_j_j(i, j, k, p) = coeff_k(p);
    _coeff_j_e(i, j, k, p) = coeff_beta(p);
    _coeff_j_sum_j(i, j, k, p) = 0.5 * (1 + _coeff_j_j(i, j, k, p));
  }

  const auto a = DrudeADEMethodCoefficient::a(epsilon_inf, constant::EPSILON_0,
                                              sum_beta, dt, sigma);
  const auto b = DrudeADEMethodCoefficient::b(epsilon_inf, constant::EPSILON_0,
                                              sum_beta, dt, sigma);

  _coeff_e_j_sum(i, j, k) = -b;

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
  auto& cexhy = calculation_param->fdtdCoefficient()->cexhy()(i, j, k);
  auto& cexhz = calculation_param->fdtdCoefficient()->cexhz()(i, j, k);

  correct_func(dy, dz, cexe, cexhy, cexhz);

  auto& ceye = calculation_param->fdtdCoefficient()->ceye()(i, j, k);
  auto& ceyhz = calculation_param->fdtdCoefficient()->ceyhz()(i, j, k);
  auto& ceyhx = calculation_param->fdtdCoefficient()->ceyhx()(i, j, k);

  correct_func(dz, dx, ceye, ceyhz, ceyhx);

  auto& ceze = calculation_param->fdtdCoefficient()->ceze()(i, j, k);
  auto& cezhx = calculation_param->fdtdCoefficient()->cezhx()(i, j, k);
  auto& cezhy = calculation_param->fdtdCoefficient()->cezhy()(i, j, k);

  correct_func(dx, dy, ceze, cezhx, cezhy);
}

}  // namespace xfdtd