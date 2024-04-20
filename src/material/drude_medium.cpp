#include <xfdtd/common/constant.h>
#include <xfdtd/material/dispersive_material.h>

#include <utility>

#include "material/dispersive_material_equation.h"
#include "updator/dispersive_material_update_method/drude_ade_method.h"

namespace xfdtd {

auto DrudeMedium::makeDrudeMedium(std::string_view name, Real eps_inf,
                                  Array1D<Real> omega_p, Array1D<Real> gamma)
    -> std::unique_ptr<DrudeMedium> {
  auto eq = std::make_shared<DrudeEqDecision>(omega_p, gamma);
  auto update_method = std::make_unique<DrudeADEMethod>(eps_inf, eq);
  return std::make_unique<DrudeMedium>(name, eps_inf, eq,
                                       std::move(update_method));
}

DrudeMedium::DrudeMedium(std::string_view name, Real epsilon_inf,
                         const std::shared_ptr<DrudeEqDecision>& eq,
                         std::unique_ptr<DrudeADEMethod> update_method,
                         ElectroMagneticProperty emp)
    : LinearDispersiveMaterial{name, epsilon_inf, eq, std::move(update_method),
                               emp},
      _drude_eq{dynamic_cast<DrudeEqDecision*>(equationPtr())} {
  if (_drude_eq == nullptr) {
    throw XFDTDLinearDispersiveMaterialEquationException(
        "DrudeMedium: invalid equation type");
  }
}

auto DrudeMedium::omegaP() const -> Array1D<Real> {
  return _drude_eq->omegaP();
}

auto DrudeMedium::gamma() const -> Array1D<Real> { return _drude_eq->gamma(); }

}  // namespace xfdtd
