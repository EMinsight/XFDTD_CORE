#include <xfdtd/common/type_define.h>
#include <xfdtd/material/dispersive_material.h>

#include <complex>
#include <memory>
#include <string_view>
#include <utility>

#include "material/dispersive_material_equation.h"
#include "updator/dispersive_material_update_method/dispersive_material_update_method.h"
#include "updator/dispersive_material_update_method/m_lorentz_ade_method.h"

namespace xfdtd {

auto LinearDispersiveMaterial::makeMLorentz(
    std::string_view name, Real epsilon_inf, Array1D<Real> a_0,
    Array1D<Real> a_1, Array1D<Real> b_0, Array1D<Real> b_1, Array1D<Real> b_2,
    ElectroMagneticProperty emp) -> std::unique_ptr<LinearDispersiveMaterial> {
  return make<LinearDispersiveMaterial::Method::ADE>(
      name, epsilon_inf,
      std::make_unique<MLorentzEqDecision>(std::move(a_0), std::move(a_1),
                                           std::move(b_0), std::move(b_1),
                                           std::move(b_2)),
      emp);
}

template <LinearDispersiveMaterial::Method method>
auto LinearDispersiveMaterial::make(
    std::string_view name, Real epsilon_inf,
    const std::shared_ptr<LinearDispersiveMaterialEquation>& equation,
    ElectroMagneticProperty emp) -> std::unique_ptr<LinearDispersiveMaterial> {
  using Method = LinearDispersiveMaterial::Method;

  if constexpr (method == Method::ADE) {
    return std::make_unique<LinearDispersiveMaterial>(
        name, epsilon_inf, equation,
        std::make_unique<MLorentzADEMethod>(
            epsilon_inf,
            std::dynamic_pointer_cast<MLorentzEqDecision>(equation)),
        emp);
  }

  // if constexpr (method == Method::DE) {
  //   static_assert(false, "Invalid method");
  // }

  return nullptr;
}

LinearDispersiveMaterial::LinearDispersiveMaterial(
    std::string_view name, Real epsilon_inf,
    std::shared_ptr<LinearDispersiveMaterialEquation> eq,
    std::unique_ptr<LinearDispersiveMaterialUpdateMethod> update_method,
    ElectroMagneticProperty emp)
    : Material{name, emp, true},
      _epsilon_inf{epsilon_inf},
      _eq{std::move(eq)},
      _update_method{std::move(update_method)} {}

LinearDispersiveMaterial::~LinearDispersiveMaterial() = default;

auto LinearDispersiveMaterial::numPoles() const -> Index {
  return _eq->numPoles();
}

auto LinearDispersiveMaterial::susceptibility(Real freq, Index p) const
    -> std::complex<Real> {
  return _eq->susceptibility(freq, p);
}

auto LinearDispersiveMaterial::relativePermittivity(
    const Array1D<Real>& freq) const -> Array1D<std::complex<Real>> {
  return xt::make_lambda_xfunction(
      [epsilon_inf = _epsilon_inf, &eq = _eq](const auto& freq) {
        auto sum{std::complex<Real>{0, 0}};
        for (Index p{0}; p < eq->numPoles(); ++p) {
          sum += eq->susceptibility(freq, p);
        }
        return epsilon_inf + sum;
      },
      freq);
}

auto MLorentzMaterial::mLorentzEquationPtr() const -> MLorentzEqDecision* {
  return dynamic_cast<MLorentzEqDecision*>(equationPtr());
}

}  // namespace xfdtd
