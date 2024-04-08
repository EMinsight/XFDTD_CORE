#include <xfdtd/material/material.h>

#include <cmath>
#include <memory>
#include <string>
#include <string_view>
#include <xtensor.hpp>

namespace xfdtd {

ElectroMagneticProperty::ElectroMagneticProperty(Real epsilon_r, Real mu_r,
                                                 Real sigma_e, Real sigma_m)
    : _epsilon_r{epsilon_r},
      _mu_r{mu_r},
      _sigma_e{sigma_e},
      _sigma_m{sigma_m},
      _epsilon{epsilon_r * EPSILON_0},
      _mu{mu_r * MU_0},
      _refract_index{std::sqrt(_epsilon_r * _mu_r)} {
  if (_sigma_e == 0) {
    _sigma_e = SIGMA_E_0;
  }
  if (_sigma_m == 0) {
    _sigma_m = SIGMA_M_0;
  }
}

ElectroMagneticProperty ElectroMagneticProperty::air() {
  return {1, 1, SIGMA_E_0, SIGMA_M_0};
}

ElectroMagneticProperty ElectroMagneticProperty::pec() {
  return {1, 1, 1e10, SIGMA_M_0};
}

ElectroMagneticProperty ElectroMagneticProperty::pmc() {
  return {1, 1, SIGMA_E_0, 1e10};
}

std::string ElectroMagneticProperty::toString() const {
  return "ElectroMagneticProperty{epsilon_r: " + std::to_string(_epsilon_r) +
         ", mu_r: " + std::to_string(_mu_r) +
         ", sigma_e: " + std::to_string(_sigma_e) +
         ", sigma_m: " + std::to_string(_sigma_m) + "}";
}

Real ElectroMagneticProperty::epsilonR() const { return _epsilon_r; }

Real ElectroMagneticProperty::muR() const { return _mu_r; }

Real ElectroMagneticProperty::sigmaE() const { return _sigma_e; }

Real ElectroMagneticProperty::sigmaM() const { return _sigma_m; }

Real ElectroMagneticProperty::epsilon() const { return _epsilon; }

Real ElectroMagneticProperty::mu() const { return _mu; }

std::complex<Real> ElectroMagneticProperty::refractIndex() const {
  return _refract_index;
}

Material::Material(std::string_view name, ElectroMagneticProperty em_property,
                   bool dispersion)
    : _name{name}, _em_property{em_property}, _dispersion{dispersion} {}

std::unique_ptr<Material> Material::createAir(std::string_view name) {
  return std::make_unique<Material>(name, ElectroMagneticProperty::air());
}

std::unique_ptr<Material> Material::createPec(std::string_view name) {
  return std::make_unique<Material>(name, ElectroMagneticProperty::pec());
}

std::unique_ptr<Material> Material::createPmc(std::string_view name) {
  return std::make_unique<Material>(name, ElectroMagneticProperty::pmc());
}

std::string Material::toString() const {
  return "Material{name: " + _name +
         ", em_property: " + _em_property.toString() +
         ", dispersion: " + std::to_string(static_cast<int>(_dispersion)) + "}";
}

std::string Material::name() const { return _name; }

ElectroMagneticProperty Material::emProperty() const { return _em_property; }

bool Material::dispersion() const { return _dispersion; }

}  // namespace xfdtd
