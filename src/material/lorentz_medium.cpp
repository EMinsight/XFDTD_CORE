#include <xfdtd/calculation_param/calculation_param.h>
#include <xfdtd/common/constant.h>
#include <xfdtd/material/dispersive_material.h>

namespace xfdtd {

template <typename T, typename F>
static auto lorentzSusceptibility(const T& omega_p, const T& eps_delta,
                                  const T& nv, const F& freq) {
  auto omega = 2 * constant::PI * freq;

  return (eps_delta * omega_p * omega_p) /
         (omega_p * omega_p + 2.0 * constant::II * nv * omega - omega * omega);
}

LorentzMedium::LorentzMedium(const std::string& name, Real eps_inf,
                             Array1D<Real> eps_static, Array1D<Real> omega_p,
                             Array1D<Real> nv)
    : LinearDispersiveMaterial{name, LinearDispersiveMaterial::Type::LORENTZ},
      _eps_inf{eps_inf},
      _eps_static{std::move(eps_static)},
      _omega_p{std::move(omega_p)},
      _nv{std::move(nv)} {
  if (_eps_static.size() != _omega_p.size() || _omega_p.size() != _nv.size()) {
    throw std::runtime_error(
        "LorentzMedium::LorentzMedium: "
        "eps_static.size() != omega_p.size() || omega_p.size() != nv.size()");
  }
}

Array1D<std::complex<Real>> LorentzMedium::relativePermittivity(
    const Array1D<Real>& freq) const {
  return xt::make_lambda_xfunction(
      [this](const auto& f) {
        std::complex<Real> sum{0, 0};
        for (std::size_t p = 0; p < numberOfPoles(); ++p) {
          sum += susceptibility(f, p);
        }
        return _eps_inf + sum;
      },
      freq);
}

std::complex<Real> LorentzMedium::susceptibility(Real freq, size_t p) const {
  if (numberOfPoles() <= p) {
    throw std::runtime_error(
        "LorentzMedium::susceptibility: "
        "p is out of range");
  }
  return lorentzSusceptibility(_omega_p(p), (_eps_static(p) - _eps_inf), _nv(p),
                               freq);
}

void LorentzMedium::calculateCoeff(const GridSpace* grid_space,
                                   const CalculationParam* calculation_param,
                                   const EMF* emf) {
  calculateCoeffForADE(calculation_param);
}

void LorentzMedium::calculateCoeffForADE(
    const CalculationParam* calculation_param) {
  auto temp_ade = [](const auto& nv_p, const auto& dt) {
    return (nv_p * dt + 1);
  };

  auto ade_alpha = [](const auto& omega_p, const auto& dt, const auto& temp) {
    return (2 - omega_p * omega_p * dt * dt) / (temp);
  };

  auto ade_xi = [](const auto& nv_p, const auto& dt, const auto& temp) {
    return (nv_p * dt - 1) / (temp);
  };

  auto ade_gamma = [eps_0 = constant::EPSILON_0](
                       const auto& omega_p, const auto& eps_delta,
                       const auto& dt, const auto& temp) {
    return (eps_0 * eps_delta * omega_p * omega_p * dt * dt) / (temp);
  };

  auto n{_omega_p.size()};

  Array1D<Real> alpha = xt::zeros<Real>({n});
  Array1D<Real> xi = xt::zeros<Real>({n});
  Array1D<Real> gamma = xt::zeros<Real>({n});
  const auto& dt = calculation_param->timeParam()->dt();
  for (std::size_t i = 0; i < n; ++i) {
    const auto& temp = temp_ade(_nv(i), dt);
    alpha(i) = ade_alpha(_omega_p(i), dt, temp);
    xi(i) = ade_xi(_nv(i), dt, temp);
    gamma(i) = ade_gamma(_omega_p(i), (_eps_static(i) - _eps_inf), dt, temp);
  }

  const auto& eps_0 = constant::EPSILON_0;
  const auto& eps_inf = _eps_inf;
  const auto& sigma_e = emProperty().sigmaE();
  const auto& sum_gamma = std::accumulate(
      gamma.begin(), gamma.end(), 0.0,
      [](Real acc, const auto& gamma_i) { return acc + gamma_i; });

  const auto& coeff_a = 2 * eps_0 * eps_inf + 0.5 * sum_gamma + sigma_e * dt;
  const auto& c1 = (0.5 * sum_gamma) / coeff_a;
  const auto& c2 = (2 * eps_0 * eps_inf - sigma_e * dt) / coeff_a;
  const auto& c3 = (2 * dt) / coeff_a;

  _coeff_for_ade = ade::LorentzCoeff{
      std::move(alpha), std::move(xi), std::move(gamma), c1, c2, c3};
}

}  // namespace xfdtd
