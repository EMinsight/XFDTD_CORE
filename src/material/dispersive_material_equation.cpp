#include "material/dispersive_material_equation.h"

#include <xfdtd/common/constant.h>
#include <xfdtd/common/index_task.h>
#include <xfdtd/common/type_define.h>
#include <xfdtd/material/dispersive_material.h>
#include <xfdtd/material/material.h>

#include <complex>
#include <utility>
#include <xtensor.hpp>

namespace xfdtd {

MLorentzEqDecision::MLorentzEqDecision(Array1D<Real> a_0, Array1D<Real> a_1,
                                       Array1D<Real> b_0, Array1D<Real> b_1,
                                       Array1D<Real> b_2)
    : _a_0{std::move(a_0)},
      _a_1{std::move(a_1)},
      _b_0{std::move(b_0)},
      _b_1{std::move(b_1)},
      _b_2{std::move(b_2)} {
  static auto anyone_not = [](auto&& key, auto&&... args) -> bool {
    return ((args != key) || ...);
  };

  if (anyone_not(a_0.size(), a_1.size(), b_0.size(), b_1.size(), b_2.size())) {
    throw XFDTDLinearDispersiveMaterialEquationException(
        "Modified Lorentz Model paramater doesn't have same size.");
  }

  _num_p = _a_0.size();
}

auto MLorentzEqDecision::susceptibility(Real freq, Index p) const
    -> std::complex<Real> {
  auto omega = 2.0 * constant::PI * freq;
  return (_a_1.at(p) * constant::II * omega + _a_0.at(p)) /
         (_b_2.at(p) * std::pow((constant::II * omega), 2) +
          _b_1.at(p) * (constant::II * omega) + _b_0.at(p));
};

DebyeEqDecision::DebyeEqDecision(Real epsilon_inf, Array1D<Real> epsilon_static,
                                 Array1D<Real> tau)
    : MLorentzEqDecision{epsilon_static - epsilon_inf,
                         xt::zeros_like(epsilon_static),
                         xt::ones_like(epsilon_static), std::move(tau),
                         xt::zeros_like(epsilon_static)},
      _epsilon_static{std::move(epsilon_static)} {}

DrudeEqDecision::DrudeEqDecision(Array1D<Real> omega_p, Array1D<Real> gamma)
    : MLorentzEqDecision{omega_p * omega_p, xt::zeros_like(omega_p),
                         xt::zeros_like(gamma), std::move(gamma),
                         xt::ones_like(omega_p)} {}

}  // namespace xfdtd