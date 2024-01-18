#ifndef _XFDTD_LIB_MATERIAL_H_
#define _XFDTD_LIB_MATERIAL_H_

#include <xfdtd/calculation_param/calculation_param.h>
#include <xfdtd/electromagnetic_field/electromagnetic_field.h>
#include <xfdtd/grid_space/grid_space.h>
#include <xfdtd/util/constant.h>

#include <complex>
#include <memory>
#include <string>

namespace xfdtd {

class ElectroMagneticProperty {
 public:
  ElectroMagneticProperty(double epsilon_r, double mu_r, double sigma_e,
                          double sigma_m);

  ElectroMagneticProperty(const ElectroMagneticProperty &) = default;

  ElectroMagneticProperty(ElectroMagneticProperty &&) noexcept = default;

  ElectroMagneticProperty &operator=(const ElectroMagneticProperty &) = default;

  ElectroMagneticProperty &operator=(ElectroMagneticProperty &&) noexcept =
      default;

  ~ElectroMagneticProperty() = default;

  static ElectroMagneticProperty air();

  static ElectroMagneticProperty pec();

  static ElectroMagneticProperty pmc();

  inline static constexpr double EPSILON_0 = constant::EPSILON_0;
  inline static constexpr double MU_0 = constant::MU_0;
  inline static constexpr double SIGMA_E_0 = 1e-20;
  inline static constexpr double SIGMA_M_0 = 1e-20;

  std::string toString() const;

  double epsilonR() const;

  double muR() const;

  double sigmaE() const;

  double sigmaM() const;

  double epsilon() const;

  double mu() const;

  std::complex<double> refractIndex() const;

 private:
  double _epsilon_r;
  double _mu_r;
  double _sigma_e;
  double _sigma_m;

  double _epsilon;
  double _mu;
  std::complex<double> _refract_index;
};

class Material {
 public:
  Material(std::string name, ElectroMagneticProperty em_property,
           bool dispersion = false);

  Material(const Material &) = default;

  Material(Material &&) noexcept = default;

  Material &operator=(const Material &) = default;

  Material &operator=(Material &&) noexcept = default;

  virtual ~Material() = default;

  static std::unique_ptr<Material> createAir(
      std::string name = "xfdtd_default_air_material");

  static std::unique_ptr<Material> createPec(
      std::string name = "xfdtd_default_pec_material");

  static std::unique_ptr<Material> createPmc(
      std::string name = "xfdtd_default_pmc_material");

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

#endif  // _XFDTD_LIB_MATERIAL_H_
