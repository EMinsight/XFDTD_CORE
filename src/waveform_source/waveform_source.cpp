#include <xfdtd/waveform_source/waveform_source.h>
#include <xfdtd/common/index_task.h>

#include <utility>

namespace xfdtd {

WaveformSource::WaveformSource(std::unique_ptr<Waveform> waveform)
    : _waveform{std::move(waveform)} {}

const std::unique_ptr<Waveform> &WaveformSource::waveform() {
  return _waveform;
}

void WaveformSource::defaultInit(
    std::shared_ptr<GridSpace> grid_space,
    std::shared_ptr<CalculationParam> calculation_param,
    std::shared_ptr<EMF> emf) {
  _grid_space = std::move(grid_space);
  _calculation_param = std::move(calculation_param);
  _emf = std::move(emf);
}

GridSpace *WaveformSource::gridSpacePtr() { return _grid_space.get(); }

CalculationParam *WaveformSource::calculationParamPtr() {
  return _calculation_param.get();
}

EMF *WaveformSource::emfPtr() { return _emf.get(); }

}  // namespace xfdtd
