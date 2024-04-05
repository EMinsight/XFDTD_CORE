#include "xfdtd/waveform_source/hard_point.h"

namespace xfdtd {

HardPoint::HardPoint(std::unique_ptr<Waveform> waveform)
    : WaveformSource(std::move(waveform)) {}

void HardPoint::init(std::shared_ptr<GridSpace> grid_space,
                     std::shared_ptr<CalculationParam> calculation_param,
                     std::shared_ptr<EMF> emf) {
  defaultInit(std::move(grid_space), std::move(calculation_param),
              std::move(emf));

  waveform()->init(calculationParamPtr()->timeParam()->eTime());
}

void HardPoint::correctMaterialSpace() {}

void HardPoint::correctUpdateCoefficient() {}

void HardPoint::updateWaveformSource() {}

// void HardPoint::correctE() {
//   // center
//   auto center{gridSpacePtr()->box().center()};
//   auto& ez{emfPtr()->ez()};
//   auto& ex{emfPtr()->ex()};
//   auto& ey{emfPtr()->ey()};

//   ez(center.i(), center.j(), center.k()) +=
//       waveform()
//           ->value()[calculationParamPtr()->timeParam()->currentTimeStep()];
//   //   ex(center.i(), center.j(), center.k()) +=
//   //       waveform()
//   //           ->value()[calculationParamPtr()->timeParam()->currentTimeStep()];

//   //   ey(center.i(), center.j(), center.k()) +=
//   //       waveform()
//   //           ->value()[calculationParamPtr()->timeParam()->currentTimeStep()];
// }

// void HardPoint::correctH() {}

std::unique_ptr<Corrector> HardPoint::generateCorrector(
    const Divider::IndexTask& task) {
  throw std::runtime_error("Not implemented");
}

}  // namespace xfdtd