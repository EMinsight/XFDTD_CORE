#ifndef _XFDTD_CORE_CAPACITOR_H_
#define _XFDTD_CORE_CAPACITOR_H_

#include "xfdtd/object/lumped_element/lumped_element.h"
namespace xfdtd {

class Capacitor : public LumpedElement {
 public:
  Capacitor(std::string name, std::unique_ptr<Cube> cube, Axis::XYZ xyz,
            Real capacitance,
            std::unique_ptr<Material> material = Material::createAir());

  Capacitor(const Capacitor&) = delete;

  Capacitor(Capacitor&&) noexcept = default;

  Capacitor& operator=(const Capacitor&) = delete;

  Capacitor& operator=(Capacitor&&) noexcept = default;

  ~Capacitor() override = default;

  std::string toString() const override;

  void init(std::shared_ptr<const GridSpace> grid_space,
            std::shared_ptr<CalculationParam> calculation_param,
            std::shared_ptr<EMF> emf) override;

  void correctUpdateCoefficient() override;

  void correctE() override;

  void correctH() override;

 private:
  Real _capacitance;
  Real _capacitance_factor;

  Array3D<Real> _da, _db, _dc;
  Array3D<Real> _beta;
};

}  // namespace xfdtd

#endif  // _XFDTD_CORE_CAPACITOR_H_
