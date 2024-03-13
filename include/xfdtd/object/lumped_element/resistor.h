#ifndef _XFDTD_CORE_RESISTOR_H_
#define _XFDTD_CORE_RESISTOR_H_

#include "xfdtd/object/lumped_element/lumped_element.h"
namespace xfdtd {

class Resistor : public LumpedElement {
 public:
  Resistor(std::string name, std::unique_ptr<Cube> cube, Axis::XYZ xyz,
           double resistance,
           std::unique_ptr<Material> material = Material::createAir());

  Resistor(const Resistor &) = delete;

  Resistor(Resistor &&) noexcept = default;

  Resistor &operator=(const Resistor &) = delete;

  Resistor &operator=(Resistor &&) noexcept = default;

  ~Resistor() override = default;

  std::string toString() const override;

  void init(std::shared_ptr<const GridSpace> grid_space,
            std::shared_ptr<CalculationParam> calculation_param,
            std::shared_ptr<EMF> emf) override;

  void correctUpdateCoefficient() override;

  void correctE() override;

  void correctH() override;

  double resistance() const;

 private:
  double _resistance;

  double _resistance_factor;
  xt::xarray<double> _da, _db, _dc;
  xt::xarray<double> _beta;
};

}  // namespace xfdtd

#endif  // _XFDTD_CORE_RESISTOR_H_
