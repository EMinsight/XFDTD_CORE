#include <xfdtd/common/constant.h>
#include <xfdtd/material/dispersive_material.h>

#include <complex>
#include <utility>

namespace xfdtd {

static auto debyeSusceptibility(const auto& eps_inf, const auto& eps_static,
                                const auto& tau, const auto& freq) {
  auto&& eps_delta = eps_static - eps_inf;
  auto&& omega = 2 * constant::PI * freq;
  return (eps_delta) / (1.0 + constant::II * omega * tau);
}

DebyeMedium::DebyeMedium(const std::string& name, double eps_inf,
                         Array1D<Real> eps_static, Array1D<Real> tau)
    : LinearDispersiveMaterial{name, Type::DEBYE},
      _eps_inf{eps_inf},
      _eps_static{std::move(eps_static)},
      _tau{std::move(tau)} {}

Array1D<std::complex<double>> DebyeMedium::relativePermittivity(
    const Array1D<Real>& freq) const {
  return xt::make_lambda_xfunction(
      [this](const auto& f) {
        std::complex<double> sum{0, 0};
        for (std::size_t p = 0; p < numberOfPoles(); ++p) {
          sum += susceptibility(f, p);
        }
        return _eps_inf + sum;
      },
      freq);
}

std::complex<double> DebyeMedium::susceptibility(double freq,
                                                 std::size_t p) const {
  if (numberOfPoles() <= p) {
    throw std::runtime_error(
        "DebyeMedium::susceptibility: "
        "p is out of range");
  }
  return debyeSusceptibility(_eps_inf, _eps_static(p), _tau(p), freq);
}

void DebyeMedium::calculateCoeff(const GridSpace* grid_space,
                                 const CalculationParam* calculation_param,
                                 const EMF* emf) {
  const auto& dt = calculation_param->timeParam()->dt();

  auto ade_k = [&dt](const auto& tau) {
    return (2 * tau - dt) / (2 * tau + dt);
  };

  auto ade_beta = [&dt, eps_0 = constant::EPSILON_0](const auto& eps_inf,
                                                     const auto& eps_static,
                                                     const auto& tau) {
    return (2 * eps_0 * (eps_static - eps_inf) * dt) / (2 * tau + dt);
  };

  auto ade_a = [](const auto& eps_inf, const auto& eps_0, const auto& sum_beta,
                  const auto& dt, const auto& sigma) {
    return (2 * eps_0 * eps_inf + sum_beta - dt * sigma) /
           (2 * eps_0 * eps_inf + sum_beta + dt * sigma);
  };

  auto ade_b = [](const auto& eps_inf, const auto& eps_0, const auto& sum_beta,
                  const auto& dt, const auto& sigma) {
    return (2 * dt) / (2 * eps_0 * eps_inf + sum_beta + dt * sigma);
  };

  const auto& k = ade_k(_tau);
  const auto& beta = ade_beta(_eps_inf, _eps_static, _tau);
  const auto& sum_beta = std::accumulate(beta.begin(), beta.end(), 0.0);
  const auto& sigma = emProperty().sigmaE();
  const auto& a = ade_a(_eps_inf, constant::EPSILON_0, sum_beta, dt, sigma);
  const auto& b = ade_b(_eps_inf, constant::EPSILON_0, sum_beta, dt, sigma);

  _coeff_for_ade = ade::DebyeCoeff{k, beta, a, b};
}

}  // namespace xfdtd
