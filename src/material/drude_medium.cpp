#include <xfdtd/common/constant.h>
#include <xfdtd/material/dispersive_material.h>
#include <xfdtd/material/dispersive_material_equation/dispersive_material_equation.h>

namespace xfdtd {

auto DrudeMedium::makeDrudeMedium(std::string_view name, Real eps_inf,
                                  Array1D<Real> omega_p, Array1D<Real> gamma)
    -> std::unique_ptr<DrudeMedium> {
  auto eq = std::make_shared<DrudeEqDecision>(omega_p, gamma);

  return std::make_unique<DrudeMedium>(name, eps_inf, eq);
}

DrudeMedium::DrudeMedium(std::string_view name, Real epsilon_inf,
                         const std::shared_ptr<DrudeEqDecision>& eq,
                         ElectroMagneticProperty emp)
    : LinearDispersiveMaterial{name, epsilon_inf, eq, emp},
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
