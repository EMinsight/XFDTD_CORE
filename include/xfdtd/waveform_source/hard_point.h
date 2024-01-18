#ifndef _XFDTD_LIB_HARD_POINT_H_
#define _XFDTD_LIB_HARD_POINT_H_

#include "xfdtd/waveform_source/waveform_source.h"
namespace xfdtd {

class HardPoint : public WaveformSource {
 public:
  explicit HardPoint(std::unique_ptr<Waveform> waveform);

  ~HardPoint() override = default;

  void init(std::shared_ptr<GridSpace> grid_space,
            std::shared_ptr<CalculationParam> calculation_param,
            std::shared_ptr<EMF> emf) override;

  void correctMaterialSpace() override;

  void correctUpdateCoefficient() override;

  void updateWaveformSourceE() override;

  void updateWaveformSourceH() override;

  void correctE() override;

  void correctH() override;
};

}  // namespace xfdtd

#endif  // _XFDTD_LIB_HARD_POINT_H_
