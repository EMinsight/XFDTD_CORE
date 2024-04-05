#ifndef _XFDTD_CORE_WAVEFORM_SOURCE_H_
#define _XFDTD_CORE_WAVEFORM_SOURCE_H_

#include <xfdtd/calculation_param/calculation_param.h>
#include <xfdtd/divider/divider.h>
#include <xfdtd/electromagnetic_field/electromagnetic_field.h>
#include <xfdtd/grid_space/grid_space.h>
#include <xfdtd/waveform/waveform.h>

#include <memory>

#include "corrector/corrector.h"

namespace xfdtd {

class WaveformSource {
 public:
  explicit WaveformSource(std::unique_ptr<Waveform> waveform);

  WaveformSource(const WaveformSource &) = delete;

  WaveformSource(WaveformSource &&) noexcept = default;

  WaveformSource &operator=(const WaveformSource &) = delete;

  WaveformSource &operator=(WaveformSource &&) noexcept = default;

  virtual ~WaveformSource() = default;

  virtual void init(std::shared_ptr<GridSpace> grid_space,
                    std::shared_ptr<CalculationParam> calculation_param,
                    std::shared_ptr<EMF> emf) = 0;

  virtual void correctMaterialSpace() = 0;

  virtual void correctUpdateCoefficient() = 0;

  virtual void initTimeDependentVariable() = 0;

  virtual void updateWaveformSource() = 0;

  const std::unique_ptr<Waveform> &waveform();

  virtual std::unique_ptr<Corrector> generateCorrector(
      const Divider::IndexTask &task) = 0;

  auto emf() const { return _emf; }

  auto calculationParam() const { return _calculation_param; }

  auto gridSpace() const { return _grid_space; }

 protected:
  virtual void defaultInit(std::shared_ptr<GridSpace> grid_space,
                           std::shared_ptr<CalculationParam> calculation_param,
                           std::shared_ptr<EMF> emf);

  GridSpace *gridSpacePtr();

  CalculationParam *calculationParamPtr();

  EMF *emfPtr();

 private:
  std::unique_ptr<Waveform> _waveform;

  std::shared_ptr<GridSpace> _grid_space;
  std::shared_ptr<CalculationParam> _calculation_param;
  std::shared_ptr<EMF> _emf;
};

}  // namespace xfdtd

#endif  // _XFDTD_CORE_WAVEFORM_SOURCE_H_
