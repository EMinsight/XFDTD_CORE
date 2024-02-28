#include "domain/domain.h"

#include "updator/updator.h"
#include "xfdtd/calculation_param/calculation_param.h"
#include "xfdtd/object/object.h"
#include "xfdtd/waveform_source/waveform_source.h"

namespace xfdtd {

// void Domain::run() {
//   while (!isCalculationDone()) {
//     updateE();

//     correctE();

//     updateH();

//     correctH();

//     communicate();

//     record();
//   }
// }

bool Domain::isCalculationDone() {
  return _calculation_param->timeParam()->endTimeStep() ==
         _calculation_param->timeParam()->currentTimeStep();
}

void Domain::updateE() { _updator->updateE(); }

void Domain::updateH() { _updator->updateH(); }

void Domain::correctE() {
  for (auto& object : _objects) {
    object->correctE();
  }

  for (auto& waveform_source : _waveform_sources) {
    waveform_source->correctE();
  }
}

void Domain::correctH() {
  for (auto& object : _objects) {
    object->correctH();
  }

  for (auto& waveform_source : _waveform_sources) {
    waveform_source->correctH();
  }
}

void Domain::communicate() {
  
}

}  // namespace xfdtd
