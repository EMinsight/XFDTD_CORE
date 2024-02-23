#ifndef _XFDTD_LIB_DISPERSIVE_MATERIAL_H_
#define _XFDTD_LIB_DISPERSIVE_MATERIAL_H_

#include <xfdtd/material/material.h>

#include <string_view>

#include "material/dispersive_solver_common.h"
#include "xfdtd/calculation_param/calculation_param.h"
#include "xfdtd/electromagnetic_field/electromagnetic_field.h"
#include "xfdtd/grid_space/grid_space.h"

namespace xfdtd {

class LinearDispersiveMaterial : public Material {
 public:
  enum class Type { UNKNOW, DEBY, DRUDE, LORENTZ };

  explicit LinearDispersiveMaterial(
      std::string_view name, Type type,
      ElectroMagneticProperty emp = ElectroMagneticProperty::air());

  LinearDispersiveMaterial(const LinearDispersiveMaterial& other) = default;

  LinearDispersiveMaterial(LinearDispersiveMaterial&&) noexcept = default;

  LinearDispersiveMaterial& operator=(const LinearDispersiveMaterial& other) =
      default;

  LinearDispersiveMaterial& operator=(
      LinearDispersiveMaterial&& other) noexcept = default;

  ~LinearDispersiveMaterial() override = default;

  Type type() const;

  virtual xt::xarray<std::complex<double>> relativePermittivity(
      const xt::xarray<double>& freq) const = 0;

  virtual void calculateCoeff(const GridSpace* grid_space,
                              const CalculationParam* calculation_param,
                              const EMF* emf) = 0;

 protected:
 private:
  Type _type;
};

class LorentzMedium : public LinearDispersiveMaterial {
 public:
  LorentzMedium(std::string_view name, double eps_inf,
                xt::xarray<double> eps_static, xt::xarray<double> omega_p,
                xt::xarray<double> nv);

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

  xt::xarray<std::complex<double>> relativePermittivity(
      const xt::xarray<double>& freq) const override;

  void calculateCoeff(const GridSpace* grid_space,
                      const CalculationParam* calculation_param,
                      const EMF* emf) override;

 private:
  double _eps_inf;
  xt::xarray<double> _eps_static, _omega_p, _nv;
  ade::LorentzCoeff _coeff_for_ade;

  void calculateCoeffForADE(const CalculationParam* calculation_param);
};

}  // namespace xfdtd

#endif  // _XFDTD_LIB_DISPERSIVE_MATERIAL_H_
