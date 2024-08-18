#include <xfdtd/common/type_define.h>
#include <xfdtd/material/dispersive_material.h>
#include <xfdtd/material/dispersive_material_equation/dispersive_material_equation.h>

#include <complex>
#include <memory>
#include <string_view>
#include <utility>

namespace xfdtd {

LinearDispersiveMaterial::LinearDispersiveMaterial(
    std::string_view name, Real epsilon_inf,
    std::shared_ptr<LinearDispersiveMaterialEquation> eq,
    ElectroMagneticProperty emp)
    : Material{name, emp, true},
      _epsilon_inf{epsilon_inf},
      _eq{std::move(eq)} {}

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

auto MLorentzMaterial::makeMLorentz(
    std::string_view name, Real epsilon_inf, Array1D<Real> a_0,
    Array1D<Real> a_1, Array1D<Real> b_0, Array1D<Real> b_1, Array1D<Real> b_2,
    ElectroMagneticProperty emp) -> std::unique_ptr<MLorentzMaterial> {
  auto eq = std::make_shared<MLorentzEqDecision>(std::move(a_0), std::move(a_1),
                                                 std::move(b_0), std::move(b_1),
                                                 std::move(b_2));
  return std::make_unique<MLorentzMaterial>(name, epsilon_inf, eq, emp);
}

auto MLorentzMaterial::mLorentzEquationPtr() const -> MLorentzEqDecision* {
  return dynamic_cast<MLorentzEqDecision*>(equationPtr());
}

}  // namespace xfdtd
