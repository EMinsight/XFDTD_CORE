#ifndef _XFDTD_CORE_DISPERSIVE_MATERIAL_H_
#define _XFDTD_CORE_DISPERSIVE_MATERIAL_H_

#include <xfdtd/calculation_param/calculation_param.h>
#include <xfdtd/common/constant.h>
#include <xfdtd/common/type_define.h>
#include <xfdtd/electromagnetic_field/electromagnetic_field.h>
#include <xfdtd/exception/exception.h>
#include <xfdtd/grid_space/grid_space.h>
#include <xfdtd/material/material.h>

#include <complex>
#include <memory>
#include <string_view>

namespace xfdtd {

class XFDTDLinearDispersiveMaterialException : public XFDTDException {
 public:
  explicit XFDTDLinearDispersiveMaterialException(const std::string& message)
      : XFDTDException(message) {}
};

class LinearDispersiveMaterialEquation;

class LinearDispersiveMaterial : public Material {
 public:
  enum class Method { ADE };

 public:
  LinearDispersiveMaterial(
      std::string_view name, Real epsilon_inf,
      std::shared_ptr<LinearDispersiveMaterialEquation> eq,
      ElectroMagneticProperty emp = ElectroMagneticProperty::air());

  LinearDispersiveMaterial(const LinearDispersiveMaterial&) = delete;

  LinearDispersiveMaterial(LinearDispersiveMaterial&&) noexcept = default;

  auto operator=(const LinearDispersiveMaterial&) -> LinearDispersiveMaterial& =
                                                         delete;

  LinearDispersiveMaterial& operator=(
      LinearDispersiveMaterial&& other) noexcept = default;

  ~LinearDispersiveMaterial() override;

  auto numPoles() const -> Index;

  auto epsilonInf() const { return _epsilon_inf; }

  auto equation() const { return _eq; }

  virtual std::complex<Real> susceptibility(Real freq, std::size_t p) const;

  virtual Array1D<std::complex<Real>> relativePermittivity(
      const Array1D<Real>& freq) const;

 protected:
  auto equationPtr() const { return _eq.get(); }

 private:
  Real _epsilon_inf;
  std::shared_ptr<LinearDispersiveMaterialEquation> _eq;
};

class MLorentzEqDecision;

class MLorentzMaterial : public LinearDispersiveMaterial {
 public:
  static auto makeMLorentz(
      std::string_view name, Real epsilon_inf, Array1D<Real> a_0,
      Array1D<Real> a_1, Array1D<Real> b_0, Array1D<Real> b_1,
      Array1D<Real> b_2,
      ElectroMagneticProperty emp = ElectroMagneticProperty::air())
      -> std::unique_ptr<MLorentzMaterial>;

 public:
  using LinearDispersiveMaterial::LinearDispersiveMaterial;

  ~MLorentzMaterial() override = default;

 protected:
  auto mLorentzEquationPtr() const -> MLorentzEqDecision*;
};

using LorentzMedium = MLorentzMaterial;

class DebyeEqDecision;

class DebyeMedium : public LinearDispersiveMaterial {
 public:
  static auto makeDebyeMedium(std::string_view name, Real epsilon_inf,
                              Array1D<Real> epsilon_static, Array1D<Real> tau)
      -> std::unique_ptr<DebyeMedium>;

 public:
  DebyeMedium(std::string_view name, Real epsilon_inf,
              const std::shared_ptr<DebyeEqDecision>& eq,
              ElectroMagneticProperty emp = ElectroMagneticProperty::air());

  DebyeMedium(DebyeMedium&&) noexcept = default;

  auto operator=(DebyeMedium&&) noexcept -> DebyeMedium& = default;

  ~DebyeMedium() override = default;

  auto tau() const -> Array1D<Real>;

  auto epsilonStatic() const -> Array1D<Real>;

 private:
  DebyeEqDecision* _debye_eq;
};

class DrudeEqDecision;

class DrudeMedium : public LinearDispersiveMaterial {
 public:
  static auto makeDrudeMedium(std::string_view name, Real eps_inf,
                              Array1D<Real> omega_p, Array1D<Real> gamma)
      -> std::unique_ptr<DrudeMedium>;

 public:
  DrudeMedium(std::string_view name, Real epsilon_inf,
              const std::shared_ptr<DrudeEqDecision>& eq,
              ElectroMagneticProperty emp = ElectroMagneticProperty::air());

  DrudeMedium(DrudeMedium&&) noexcept = default;

  auto operator=(DrudeMedium&&) noexcept -> DrudeMedium& = default;

  ~DrudeMedium() override = default;

  auto omegaP() const -> Array1D<Real>;

  auto gamma() const -> Array1D<Real>;

 private:
  DrudeEqDecision* _drude_eq;
};

}  // namespace xfdtd

#endif  // _XFDTD_CORE_DISPERSIVE_MATERIAL_H_
