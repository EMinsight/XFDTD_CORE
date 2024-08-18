#include <xfdtd/common/constant.h>
#include <xfdtd/common/type_define.h>
#include <xfdtd/material/dispersive_material.h>
#include <xfdtd/material/dispersive_material_equation/dispersive_material_equation.h>

#include <memory>
#include <utility>

namespace xfdtd {

auto DebyeMedium::makeDebyeMedium(std::string_view name, Real epsilon_inf,
                                  Array1D<Real> epsilon_static,
                                  Array1D<Real> tau)

    -> std::unique_ptr<DebyeMedium> {
  auto debye_eq = std::make_shared<DebyeEqDecision>(
      epsilon_inf, std::move(epsilon_static), std::move(tau));

  return std::make_unique<DebyeMedium>(name, epsilon_inf, debye_eq);
}

DebyeMedium::DebyeMedium(std::string_view name, Real epsilon_inf,
                         const std::shared_ptr<DebyeEqDecision>& eq,
                         ElectroMagneticProperty emp)
    : LinearDispersiveMaterial{name, epsilon_inf, eq, emp},
      _debye_eq{dynamic_cast<DebyeEqDecision*>(equationPtr())} {
  if (_debye_eq == nullptr) {
    throw XFDTDLinearDispersiveMaterialEquationException(
        "DebyeMedium: invalid equation type");
  }
}

auto DebyeMedium::tau() const -> Array1D<Real> { return _debye_eq->tau(); }

auto DebyeMedium::epsilonStatic() const -> Array1D<Real> {
  return _debye_eq->epsilonStatic();
}

}  // namespace xfdtd
