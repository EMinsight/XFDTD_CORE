#ifndef _XFDTD_CORE_INDUCTOR_H_
#define _XFDTD_CORE_INDUCTOR_H_

#include "xfdtd/object/lumped_element/lumped_element.h"

namespace xfdtd {

class Inductor : public LumpedElement {
 public:
  Inductor(std::string name, std::unique_ptr<Cube> cube, Axis::XYZ xyz,
           Real inductance,
           std::unique_ptr<Material> material = Material::createAir());

  Inductor(const Inductor&) = delete;

  Inductor(Inductor&&) noexcept = default;

  Inductor& operator=(const Inductor&) = delete;

  Inductor& operator=(Inductor&&) noexcept = default;

  ~Inductor() override = default;

  std::string toString() const override;

  void init(std::shared_ptr<const GridSpace> grid_space,
            std::shared_ptr<CalculationParam> calculation_param,
            std::shared_ptr<EMF> emf) override;

  void correctUpdateCoefficient() override;

  void correctE() override;

  void correctH() override;

  std::unique_ptr<Corrector> generateCorrector(
      const Task<std::size_t>& task) override;

 private:
  Real _inductance;
  Real _inductance_factor;

  Array3D<Real> _da, _db, _dc;
  Array3D<Real> _beta;
  Array3D<Real> _cecjc, _cjcec;
  Array3D<Real> _j;
};

}  // namespace xfdtd

#endif  // _XFDTD_CORE_INDUCTOR_H_
