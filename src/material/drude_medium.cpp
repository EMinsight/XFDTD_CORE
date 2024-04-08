#include <complex>
#include <numeric>

#include "xfdtd/material/dispersive_material.h"
#include "xfdtd/common/constant.h"

namespace xfdtd {

template <typename T, typename F>
static auto drudeSusceptibility(const T& omega_p, const T& gamma,
                                const F& freq) {
  using namespace std::complex_literals;
  auto omega = 2 * constant::PI * freq;

  return -(omega_p * omega_p) / (omega * omega - 1i * gamma * omega);
}

DrudeMedium::DrudeMedium(std::string_view name, double eps_inf,
                         xt::xarray<double> omega_p, xt::xarray<double> gamma)
    : LinearDispersiveMaterial{name, Type::DRUDE},
      _eps_inf{eps_inf},
      _omega_p{std::move(omega_p)},
      _gamma{std::move(gamma)} {}

xt::xarray<std::complex<double>> DrudeMedium::relativePermittivity(
    const xt::xarray<double>& freq) const {
  return xt::make_lambda_xfunction(
      [eps_inf = _eps_inf, this](const auto& f) {
        std::complex<double> sum{0, 0};
        for (std::size_t p = 0; p < numberOfPoles(); ++p) {
          sum += susceptibility(f, p);
        }
        return eps_inf + sum;
      },
      freq);
}

std::complex<double> DrudeMedium::susceptibility(double freq,
                                                 std::size_t p) const {
  if (numberOfPoles() <= p) {
    throw std::runtime_error(
        "DrudeMedium::susceptibility: "
        "p is out of range");
  }

  return drudeSusceptibility(_omega_p(p), _gamma(p), freq);
}

void DrudeMedium::calculateCoeff(const GridSpace* grid_space,
                                 const CalculationParam* calculation_param,
                                 const EMF* emf) {
  const auto& dt = calculation_param->timeParam()->dt();

  auto ade_k = [&dt](const auto& gamma) {
    return (1 - gamma * dt / 2) / (1 + gamma * dt / 2);
  };

  auto ade_beta = [&dt, eps_0 = constant::EPSILON_0](const auto& omega_p,
                                                     const auto& gamma) {
    return (eps_0 * omega_p * omega_p * dt * dt / 2) / (1 + gamma * dt / 2);
  };

  auto ade_a = [](const auto& eps_inf, const auto& eps_0, const auto& sum_beta,
                  const auto& dt, const auto& sigma) {
    return (2 * eps_0 * eps_inf - dt * sum_beta - dt * sigma) /
           (2 * eps_0 * eps_inf + dt * sum_beta + dt * sigma);
  };

  auto ade_b = [](const auto& eps_inf, const auto& eps_0, const auto& sum_beta,
                  const auto& dt, const auto& sigma) {
    return (2 * dt) / (2 * eps_0 * eps_inf + dt * sum_beta + dt * sigma);
  };

  const auto& k = ade_k(_gamma);
  const auto& beta = ade_beta(_omega_p, _gamma);

  const auto& sigma_e = emProperty().sigmaE();
  const auto& eps_0 = constant::EPSILON_0;
  const auto& eps_inf = _eps_inf;
  const auto& sum_beta = std::accumulate(beta.begin(), beta.end(), 0.0);

  const auto& a = ade_a(eps_inf, eps_0, sum_beta, dt, sigma_e);
  const auto& b = ade_b(eps_inf, eps_0, sum_beta, dt, sigma_e);

  _coeff_for_ade = ade::DrudeCoeff{k, beta, a, b};
}

}  // namespace xfdtd
