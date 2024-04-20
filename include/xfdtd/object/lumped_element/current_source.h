#ifndef _XFDTD_CORE_CURRENT_SOURCE_H_
#define _XFDTD_CORE_CURRENT_SOURCE_H_

#include "xfdtd/object/lumped_element/lumped_element.h"
#include "xfdtd/waveform/waveform.h"

namespace xfdtd {

class CurrentSource : public LumpedElement {
 public:
  CurrentSource(std::string name, std::unique_ptr<Cube> cube,
                Axis::Direction direction, Real resistance,
                std::unique_ptr<Waveform> waveform,
                std::unique_ptr<Material> material = Material::createAir());

  CurrentSource(const CurrentSource &) = delete;

  CurrentSource(CurrentSource &&) noexcept = default;

  CurrentSource &operator=(const CurrentSource &) = delete;

  CurrentSource &operator=(CurrentSource &&) noexcept = default;

  ~CurrentSource() override = default;

  std::string toString() const override;

  void init(std::shared_ptr<const GridSpace> grid_space,
            std::shared_ptr<CalculationParam> calculation_param,
            std::shared_ptr<EMF> emf) override;

  void correctUpdateCoefficient() override;

  void initTimeDependentVariable() override;

  void correctE() override;

  void correctH() override;

  std::unique_ptr<Corrector> generateCorrector(
      const Task<std::size_t> &task) override;

  Axis::Direction direction() const;

  Real resistance() const;

  const std::unique_ptr<Waveform> &waveform() const;

 private:
  Axis::Direction _direction;
  Real _resistance;
  std::unique_ptr<Waveform> _waveform;


  Real _resistance_factor;
  Real _current_amplitude_factor;
  Array3D<Real> _da, _db, _dc;
  Array3D<Real> _alpha, _beta;
  Array3D<Real> _coff_i;
};

}  // namespace xfdtd

#endif  // _XFDTD_CORE_CURRENT_SOURCE_H_
