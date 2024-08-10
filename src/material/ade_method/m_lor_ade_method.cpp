#include <xfdtd/common/constant.h>
#include <xfdtd/material/ade_method/m_lor_ade_method.h>
#include <xfdtd/material/dispersive_material_equation/dispersive_material_equation.h>

#include <xtensor.hpp>

namespace xfdtd {

class MLorentzADEMethodCoefficient {
 public:
  static auto a(const Array1D<Real>& a0, const Array1D<Real>& a1,
                Real dt) -> Array1D<Real> {
    return ((a1 * constant::EPSILON_0 / (dt * dt)) +
            (0.5 * a0 * constant::EPSILON_0 / dt));
  }

  static auto b(const Array1D<Real>& a1, Real dt) -> Array1D<Real> {
    return (-2 * a1 * constant::EPSILON_0) / (dt * dt);
  }

  static auto c(const Array1D<Real>& a0, const Array1D<Real>& a1,
                Real dt) -> Array1D<Real> {
    return ((a1 * constant::EPSILON_0 / (dt * dt)) -
            (0.5 * a0 * constant::EPSILON_0 / dt));
  }

  static auto d(const Array1D<Real>& b0, const Array1D<Real>& b2,
                Real dt) -> Array1D<Real> {
    return ((2 * b2 / (dt * dt)) - b0);
  }

  static auto e(const Array1D<Real>& b1, const Array1D<Real>& b2,
                Real dt) -> Array1D<Real> {
    return ((-b2 / (dt * dt)) + (b1 / (2 * dt)));
  }

  static auto f(const Array1D<Real>& b1, const Array1D<Real>& b2,
                Real dt) -> Array1D<Real> {
    return ((b2 / (dt * dt)) + (b1 / (2 * dt)));
  }
};

MLorentzADEMethodStorage::MLorentzADEMethodStorage(Index num_pole, Index nx,
                                                   Index ny, Index nz)
    : ADEMethodStorage{num_pole} {
  _coeff_j_j = xt::zeros<Real>({nx, ny, nz, num_pole});
  _coeff_j_j_p = xt::zeros<Real>({nx, ny, nz, num_pole});
  _coeff_j_e_n = xt::zeros<Real>({nx, ny, nz, num_pole});
  _coeff_j_e = xt::zeros<Real>({nx, ny, nz, num_pole});
  _coeff_j_e_p = xt::zeros<Real>({nx, ny, nz, num_pole});
  _coeff_j_sum_j = xt::zeros<Real>({nx, ny, nz, num_pole});
  _coeff_j_sum_j_p = xt::zeros<Real>({nx, ny, nz, num_pole});
  _coeff_e_j_sum = xt::zeros<Real>({nx, ny, nz});
  _coeff_e_e_p = xt::zeros<Real>({nx, ny, nz});
  _ex_prev = xt::zeros<Real>({nx, ny, nz});
  _ey_prev = xt::zeros<Real>({nx, ny, nz});
  _ez_prev = xt::zeros<Real>({nx, ny, nz});
  _jx_arr = xt::zeros<Real>({nx, ny, nz, num_pole});
  _jy_arr = xt::zeros<Real>({nx, ny, nz, num_pole});
  _jz_arr = xt::zeros<Real>({nx, ny, nz, num_pole});
  _jx_prev_arr = xt::zeros<Real>({nx, ny, nz, num_pole});
  _jy_prev_arr = xt::zeros<Real>({nx, ny, nz, num_pole});
  _jz_prev_arr = xt::zeros<Real>({nx, ny, nz, num_pole});
}

auto MLorentzADEMethodStorage::correctCoeff(
    Index i, Index j, Index k,
    const LinearDispersiveMaterial& linear_dispersive_material,
    const std::shared_ptr<const GridSpace>& grid_space,
    const std::shared_ptr<CalculationParam>& calculation_param) -> void {
  const auto equation = linear_dispersive_material.equation();
  const auto m_lor_eq = std::dynamic_pointer_cast<MLorentzEqDecision>(equation);
  const auto dt = calculation_param->timeParam()->dt();

  if (_num_pole < linear_dispersive_material.numPoles()) {
    std::stringstream ss;
    ss << "The number of poles is not enough. The number of poles is "
       << _num_pole << " but the number of poles in the material is "
       << linear_dispersive_material.numPoles();
    throw XFDTDLinearDispersiveMaterialException{ss.str()};
  }

  if (m_lor_eq == nullptr) {
    throw XFDTDLinearDispersiveMaterialException{
        "The equation is not MLorentzEqDecision"};
  }

  auto a = MLorentzADEMethodCoefficient::a(m_lor_eq->a0(), m_lor_eq->a1(), dt);
  auto b = MLorentzADEMethodCoefficient::b(m_lor_eq->a1(), dt);
  auto c = MLorentzADEMethodCoefficient::c(m_lor_eq->a0(), m_lor_eq->a1(), dt);
  auto d = MLorentzADEMethodCoefficient::d(m_lor_eq->b0(), m_lor_eq->b2(), dt);
  auto e = MLorentzADEMethodCoefficient::e(m_lor_eq->b1(), m_lor_eq->b2(), dt);
  auto f = MLorentzADEMethodCoefficient::f(m_lor_eq->b1(), m_lor_eq->b2(), dt);

  for (Index p{0}; p < linear_dispersive_material.numPoles(); ++p) {
    _coeff_j_e_n(i, j, k, p) = a(p) / f(p);
    _coeff_j_e(i, j, k, p) = b(p) / f(p);
    _coeff_j_e_p(i, j, k, p) = c(p) / f(p);
    _coeff_j_j(i, j, k, p) = d(p) / f(p);
    _coeff_j_j_p(i, j, k, p) = e(p) / f(p);
    _coeff_j_sum_j(i, j, k, p) = 0.5 * (1 + d(p) / f(p));
    _coeff_j_sum_j_p(i, j, k, p) = 0.5 * (e(p) / f(p));
  }

  auto sum_a_f = Real{0};
  for (Index p{0}; p < linear_dispersive_material.numPoles(); ++p) {
    sum_a_f += a(p) / f(p);
  }
  auto sum_b_f = Real{0};
  for (Index p{0}; p < linear_dispersive_material.numPoles(); ++p) {
    sum_b_f += b(p) / f(p);
  }
  Real sum_c_f = 0;
  for (Index p{0}; p < linear_dispersive_material.numPoles(); ++p) {
    sum_c_f += c(p) / f(p);
  }

  auto epsilon_inf = linear_dispersive_material.epsilonInf();
  auto sigma_e = linear_dispersive_material.emProperty().sigmaE();

  Real temp =
      (epsilon_inf * constant::EPSILON_0) / dt + 0.5 * sum_a_f + 0.5 * sigma_e;

  _coeff_e_j_sum(i, j, k) = -1 / temp;
  _coeff_e_e_p(i, j, k) = -0.5 * sum_c_f / temp;

  auto correct_func = [epsilon_inf, sum_a_f, sum_b_f, temp](
                          const auto& dt, const auto& da, const auto& db,
                          const auto& sigma_e, auto& cece, auto& cecha,
                          auto& cechb) {
    const auto epsilon_0 = xfdtd::constant::EPSILON_0;
    cece =
        ((epsilon_inf * epsilon_0) / dt - 0.5 * sum_b_f - 0.5 * sigma_e) / temp;
    cecha = -1 / (temp * db);
    cechb = 1 / (temp * da);
  };

  const auto dx = grid_space->hSizeX()(i);
  const auto dy = grid_space->hSizeY()(j);
  const auto dz = grid_space->hSizeZ()(k);

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

}  // namespace xfdtd