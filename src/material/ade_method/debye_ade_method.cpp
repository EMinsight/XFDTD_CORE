#include <xfdtd/material/ade_method/debye_ade_method.h>
#include <xfdtd/material/dispersive_material.h>
#include <xfdtd/material/dispersive_material_equation/dispersive_material_equation.h>

#include <sstream>

namespace xfdtd {

class DebyeADEMethodCoefficient {
 public:
  static auto k(const Array1D<Real>& tau, Real dt) -> Array1D<Real> {
    return (2 * tau - dt) / (2 * tau + dt);
  }

  static auto beta(const Array1D<Real>& delta_epsilon, const Array1D<Real>& tau,
                   Real dt) -> Array1D<Real> {
    return (2 * constant::EPSILON_0 * delta_epsilon * dt) / (2 * tau + dt);
  }

  static auto a(const auto& epsilon_inf, const auto& epsilon_0,
                const auto& sum_beta, const auto& dt, const auto& sigma) {
    return (2 * epsilon_0 * epsilon_inf + sum_beta - dt * sigma) /
           (2 * epsilon_0 * epsilon_inf + sum_beta + dt * sigma);
  }

  static auto b(const auto& epsilon_inf, const auto& epsilon_0,
                const auto& sum_beta, const auto& dt, const auto& sigma) {
    return (2 * dt) / (2 * epsilon_0 * epsilon_inf + sum_beta + dt * sigma);
  }
};

DebyeADEMethodStorage::DebyeADEMethodStorage(Index num_pole, Index nx, Index ny,
                                             Index nz)
    : ADEMethodStorage{num_pole} {
  _coeff_j_j = xt::zeros<Real>({nx, ny, nz, num_pole});
  _coeff_j_e = xt::zeros<Real>({nx, ny, nz, num_pole});
  _coeff_j_sum_j = xt::zeros<Real>({nx, ny, nz, num_pole});
  _coeff_e_j_sum = xt::zeros<Real>({nx, ny, nz});
  _ex_prev = xt::zeros<Real>({nx, ny, nz});
  _ey_prev = xt::zeros<Real>({nx, ny, nz});
  _ez_prev = xt::zeros<Real>({nx, ny, nz});
  _jx_arr = xt::zeros<Real>({nx, ny, nz, num_pole});
  _jy_arr = xt::zeros<Real>({nx, ny, nz, num_pole});
  _jz_arr = xt::zeros<Real>({nx, ny, nz, num_pole});
}

auto DebyeADEMethodStorage::correctCoeff(
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

  auto debye_eq = std::dynamic_pointer_cast<DebyeEqDecision>(
      linear_dispersive_material.equation());
  if (!debye_eq) {
    throw XFDTDLinearDispersiveMaterialException{
        "The equation is not Debye equation"};
  }

  auto dt = calculation_param->timeParam()->dt();

  auto delta_epsilon = debye_eq->deltaEpsilon();
  auto tau = debye_eq->tau();

  auto coeff_k = DebyeADEMethodCoefficient::k(tau, dt);
  auto coeff_beta = DebyeADEMethodCoefficient::beta(delta_epsilon, tau, dt);
  auto sum_beta =
      std::accumulate(coeff_beta.begin(), coeff_beta.end(), Real{0.0});

  for (Index p{0}; p < _num_pole; ++p) {
    _coeff_j_j(i, j, k, p) = coeff_k[p];
    _coeff_j_e(i, j, k, p) = coeff_beta[p] / dt;
    _coeff_j_sum_j(i, j, k, p) = 0.5 * (1 + coeff_k[p]);
  }

  auto epsilon_inf = linear_dispersive_material.epsilonInf();
  auto sigma_e = linear_dispersive_material.emProperty().sigmaE();

  const auto temp_a = DebyeADEMethodCoefficient::a(
      epsilon_inf, constant::EPSILON_0, sum_beta, dt, sigma_e);
  const auto temp_b = DebyeADEMethodCoefficient::b(
      epsilon_inf, constant::EPSILON_0, sum_beta, dt, sigma_e);

  _coeff_e_j_sum(i, j, k) = -temp_b;

  auto correct_func = [a = temp_a, b = temp_b](const auto da, const auto db,
                                               auto& cece, auto& cecha,
                                               auto& cechb) {
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