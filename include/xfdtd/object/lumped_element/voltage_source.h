#ifndef _XFDTD_CORE_VOLTAGE_SOURCE_H_
#define _XFDTD_CORE_VOLTAGE_SOURCE_H_

#include <memory>

#include "xfdtd/coordinate_system/coordinate_system.h"
#include "xfdtd/grid_space/grid_space.h"
#include "xfdtd/object/lumped_element/lumped_element.h"
#include "xfdtd/waveform/waveform.h"
namespace xfdtd {

class VoltageSource : public LumpedElement {
 public:
  VoltageSource(std::string name, std::unique_ptr<Cube> cube,
                Axis::Direction direction, double resistance,
                std::unique_ptr<Waveform> waveform,
                std::unique_ptr<Material> material = Material::createAir());

  VoltageSource(const VoltageSource &) = delete;

  VoltageSource(VoltageSource &&) noexcept = default;

  VoltageSource &operator=(const VoltageSource &) = delete;

  VoltageSource &operator=(VoltageSource &&) noexcept = default;

  ~VoltageSource() override = default;

  std::string toString() const override;

  void init(std::shared_ptr<const GridSpace> grid_space,
            std::shared_ptr<CalculationParam> calculation_param,
            std::shared_ptr<EMF> emf) override;

  void correctUpdateCoefficient() override;

  void initTimeDependentVariable() override;

  void correctE() override;

  void correctH() override;

  std::unique_ptr<Corrector> generateCorrector(
      const Divider::Task<std::size_t> &task) override;

  Axis::Direction direction() const;

  double resistance() const;

  const std::unique_ptr<Waveform> &waveform() const;

 private:
  Axis::Direction _direction;
  double _resistance;
  std::unique_ptr<Waveform> _waveform;

  double _resistance_factor;
  double _voltage_amplitude_factor;
  xt::xarray<double> _da, _db, _dc;
  xt::xarray<double> _alpha, _beta;
  xt::xarray<double> _coff_v;
};

}  // namespace xfdtd

#endif  // _XFDTD_CORE_VOLTAGE_SOURCE_H_
