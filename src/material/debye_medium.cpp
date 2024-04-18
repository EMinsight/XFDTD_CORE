#include <xfdtd/common/constant.h>
#include <xfdtd/common/type_define.h>
#include <xfdtd/material/dispersive_material.h>

#include <memory>
#include <utility>

#include "material/dispersive_material_equation.h"
#include "updator/dispersive_material_update_method/debye_ade_method.h"

namespace xfdtd {

auto DebyeMedium::makeDebyeMedium(std::string_view name, Real epsilon_inf,
                                  Array1D<Real> epsilon_static,
                                  Array1D<Real> tau)

    -> std::unique_ptr<DebyeMedium> {
  auto debye_eq = std::make_shared<DebyeEqDecision>(
      epsilon_inf, std::move(epsilon_static), std::move(tau));

  auto update_method = std::make_unique<DebyeADEMethod>(epsilon_inf, debye_eq);
  return std::make_unique<DebyeMedium>(name, epsilon_inf, debye_eq,
                                       std::move(update_method));
}

DebyeMedium::DebyeMedium(std::string_view name, Real epsilon_inf,
                         const std::shared_ptr<DebyeEqDecision>& eq,
                         std::unique_ptr<DebyeADEMethod> update_method,
                         ElectroMagneticProperty emp)
    : LinearDispersiveMaterial{name, epsilon_inf, eq, std::move(update_method),
                               emp},
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
