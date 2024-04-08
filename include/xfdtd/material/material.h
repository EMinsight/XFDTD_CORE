#ifndef _XFDTD_CORE_MATERIAL_H_
#define _XFDTD_CORE_MATERIAL_H_

#include <xfdtd/calculation_param/calculation_param.h>
#include <xfdtd/common/constant.h>
#include <xfdtd/common/type_define.h>
#include <xfdtd/electromagnetic_field/electromagnetic_field.h>
#include <xfdtd/grid_space/grid_space.h>

#include <complex>
#include <memory>
#include <string_view>

namespace xfdtd {

class ElectroMagneticProperty {
 public:
  ElectroMagneticProperty(Real epsilon_r, Real mu_r, Real sigma_e,
                          Real sigma_m);

  ElectroMagneticProperty(const ElectroMagneticProperty &) = default;

  ElectroMagneticProperty(ElectroMagneticProperty &&) noexcept = default;

  ElectroMagneticProperty &operator=(const ElectroMagneticProperty &) = default;

  ElectroMagneticProperty &operator=(ElectroMagneticProperty &&) noexcept =
      default;

  ~ElectroMagneticProperty() = default;

  static ElectroMagneticProperty air();

  static ElectroMagneticProperty pec();

  static ElectroMagneticProperty pmc();

  inline static constexpr Real EPSILON_0 = constant::EPSILON_0;
  inline static constexpr Real MU_0 = constant::MU_0;
  inline static constexpr Real SIGMA_E_0 = 1e-20;
  inline static constexpr Real SIGMA_M_0 = 1e-20;

  std::string toString() const;

  Real epsilonR() const;

  Real muR() const;

  Real sigmaE() const;

  Real sigmaM() const;

  Real epsilon() const;

  Real mu() const;

  std::complex<Real> refractIndex() const;

 private:
  Real _epsilon_r;
  Real _mu_r;
  Real _sigma_e;
  Real _sigma_m;

  Real _epsilon;
  Real _mu;
  std::complex<Real> _refract_index;
};

class Material {
 public:
  Material(std::string_view name, ElectroMagneticProperty em_property,
           bool dispersion = false);

  Material(const Material &) = default;

  Material(Material &&) noexcept = default;

  Material &operator=(const Material &) = default;

  Material &operator=(Material &&) noexcept = default;

  virtual ~Material() = default;

  static std::unique_ptr<Material> createAir(
      std::string_view name = "xfdtd_default_air_material");

  static std::unique_ptr<Material> createPec(
      std::string_view name = "xfdtd_default_pec_material");

  static std::unique_ptr<Material> createPmc(
      std::string_view name = "xfdtd_default_pmc_material");

  virtual std::string toString() const;

  std::string name() const;

  ElectroMagneticProperty emProperty() const;

  bool dispersion() const;

 private:
  std::string _name;
  ElectroMagneticProperty _em_property;
  bool _dispersion;
};

}  // namespace xfdtd

#endif  // _XFDTD_CORE_MATERIAL_H_
