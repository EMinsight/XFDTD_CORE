#ifndef _XFDTD_CORE_DISPERSIVE_MATERIAL_H_
#define _XFDTD_CORE_DISPERSIVE_MATERIAL_H_

#include <xfdtd/calculation_param/calculation_param.h>
#include <xfdtd/common/type_define.h>
#include <xfdtd/electromagnetic_field/electromagnetic_field.h>
#include <xfdtd/grid_space/grid_space.h>
#include <xfdtd/material/material.h>

namespace xfdtd {

namespace ade {
struct LorentzCoeff {
  Array1D<Real> _alpha;
  Array1D<Real> _xi;
  Array1D<Real> _gamma;

  Real _c1;
  Real _c2;
  Real _c3;
};

struct DrudeCoeff {
  Array1D<Real> _k;
  Array1D<Real> _beta;

  Real _a;
  Real _b;
};

struct DebyCoeff {
  Array1D<Real> _k;
  Array1D<Real> _beta;

  Real _a;
  Real _b;
};

}  // namespace ade

class LinearDispersiveMaterial : public Material {
 public:
  enum class Type { UNKNOW, DEBYE, DRUDE, LORENTZ };

  explicit LinearDispersiveMaterial(
      const std::string& name, Type type,
      ElectroMagneticProperty emp = ElectroMagneticProperty::air());

  LinearDispersiveMaterial(const LinearDispersiveMaterial& other) = default;

  LinearDispersiveMaterial(LinearDispersiveMaterial&&) noexcept = default;

  LinearDispersiveMaterial& operator=(const LinearDispersiveMaterial& other) =
      default;

  LinearDispersiveMaterial& operator=(
      LinearDispersiveMaterial&& other) noexcept = default;

  ~LinearDispersiveMaterial() override = default;

  Type type() const;

  virtual Array1D<std::complex<Real>> relativePermittivity(
      const Array1D<Real>& freq) const = 0;

  virtual void calculateCoeff(const GridSpace* grid_space,
                              const CalculationParam* calculation_param,
                              const EMF* emf) = 0;

  virtual std::complex<Real> susceptibility(Real freq, std::size_t p) const = 0;

 private:
  Type _type;
};

class LorentzMedium : public LinearDispersiveMaterial {
 public:
  LorentzMedium(const std::string& name, Real eps_inf,
                Array1D<Real> eps_static, Array1D<Real> omega_p,
                Array1D<Real> nv);

  LorentzMedium(const LorentzMedium& other) = default;

  LorentzMedium(LorentzMedium&& other) noexcept = default;

  LorentzMedium& operator=(const LorentzMedium& other) = default;

  LorentzMedium& operator=(LorentzMedium&& other) noexcept = default;

  ~LorentzMedium() override = default;

  auto numberOfPoles() const { return _omega_p.size(); }

  auto epsStatic() const { return _eps_static; }

  auto epsInf() const { return _eps_inf; }

  auto omegaP() const { return _omega_p; }

  auto nv() const { return _nv; }

  const auto& coeffForADE() const { return _coeff_for_ade; }

  Array1D<std::complex<Real>> relativePermittivity(
      const Array1D<Real>& freq) const override;

  std::complex<Real> susceptibility(Real freq, std::size_t p) const override;

  void calculateCoeff(const GridSpace* grid_space,
                      const CalculationParam* calculation_param,
                      const EMF* emf) override;

 private:
  Real _eps_inf;
  Array1D<Real> _eps_static, _omega_p, _nv;
  ade::LorentzCoeff _coeff_for_ade;

  void calculateCoeffForADE(const CalculationParam* calculation_param);
};

class DrudeMedium : public LinearDispersiveMaterial {
 public:
  DrudeMedium(const std::string& name, Real eps_inf, Array1D<Real> omega_p,
              Array1D<Real> gamma);

  ~DrudeMedium() override = default;

  auto numberOfPoles() const { return _omega_p.size(); }

  auto epsInf() const { return _eps_inf; }

  auto omegaP() const { return _omega_p; }

  auto gamma() const { return _gamma; }

  const auto& coeffForADE() const { return _coeff_for_ade; }

  Array1D<std::complex<Real>> relativePermittivity(
      const Array1D<Real>& freq) const override;

  std::complex<Real> susceptibility(Real freq, std::size_t p) const override;

  void calculateCoeff(const GridSpace* grid_space,
                      const CalculationParam* calculation_param,
                      const EMF* emf) override;

 private:
  Real _eps_inf;
  Array1D<Real> _omega_p, _gamma;
  ade::DrudeCoeff _coeff_for_ade;

  void calculateCoeffForADE(const CalculationParam* calculation_param);
};

class DebyeMedium : public LinearDispersiveMaterial {
 public:
  DebyeMedium(const std::string& name, Real eps_inf,
              Array1D<Real> eps_static, Array1D<Real> tau);

  ~DebyeMedium() override = default;

  auto numberOfPoles() const { return _tau.size(); }

  auto tau() const { return _tau; }

  auto epsInf() const { return _eps_inf; }

  const auto& epsStatic() const { return _eps_static; }

  const auto& coeffForADE() const { return _coeff_for_ade; }

  Array1D<std::complex<Real>> relativePermittivity(
      const Array1D<Real>& freq) const override;

  std::complex<Real> susceptibility(Real freq, std::size_t p) const override;

  void calculateCoeff(const GridSpace* grid_space,
                      const CalculationParam* calculation_param,
                      const EMF* emf) override;

 private:
  Real _eps_inf;
  Array1D<Real> _eps_static, _tau;
  ade::DebyCoeff _coeff_for_ade;
};

}  // namespace xfdtd

#endif  // _XFDTD_CORE_DISPERSIVE_MATERIAL_H_
