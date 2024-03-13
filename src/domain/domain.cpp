#include "domain/domain.h"

#include <mutex>
#include <ostream>
#include <sstream>
#include <xtensor/xbuilder.hpp>

#include "updator/updator.h"

namespace xfdtd {

static std::mutex cout_mutex;

void Domain::run() {
  while (!isCalculationDone()) {
    updateE();

    correctE();

    updateH();

    correctH();

    record();

    nextStep();
  }
}

bool Domain::isCalculationDone() const {
  return _calculation_param->timeParam()->endTimeStep() <=
         _calculation_param->timeParam()->currentTimeStep();
}

void Domain::updateE() {
  communicateForH();
  _updator->updateE();
  if (isMaster()) {
    for (auto&& ws : _waveform_sources) {
      ws->updateWaveformSourceE();
    }
  }
  synchronize();
}

void Domain::updateH() {
  communicateForE();
  _updator->updateH();
  if (isMaster()) {
    for (auto&& ws : _waveform_sources) {
      ws->updateWaveformSourceH();
    }
  }
  synchronize();
}

void Domain::correctE() {
  for (auto&& c : _correctors) {
    c->correctE();
  }
  synchronize();
}

void Domain::correctH() {
  for (auto&& c : _correctors) {
    c->correctH();
  }
  synchronize();
}

void Domain::synchronize() { _barrier.arrive_and_wait(); }

void Domain::communicate() { synchronize(); }

void Domain::record() {
  for (auto&& m : _monitors) {
    m->update();
  }

  for (auto&& n : _nfffts) {
    n->update();
  }

  synchronize();
}

void Domain::nextStep() {
  if (isMaster()) {
    std::stringstream ss;
    ss << "\r"
       << "Progress: " << _calculation_param->timeParam()->currentTimeStep() + 1
       << "/" << _calculation_param->timeParam()->endTimeStep();
    std::cout << ss.str() << std::flush;
    _calculation_param->timeParam()->nextStep();
  }
  synchronize();
}

void Domain::communicateForH() {}

void Domain::communicateForE() {}

}  // namespace xfdtd
