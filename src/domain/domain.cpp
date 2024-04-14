#include "domain/domain.h"

#include <xfdtd/parallel/mpi_support.h>
#include <xfdtd/util/fdtd_basic.h>

#include <cstdlib>
#include <ios>
#include <mutex>
#include <ostream>
#include <sstream>
#include <string>
#include <xtensor.hpp>
#include <xtensor/xbuilder.hpp>

#include "updator/updator.h"

namespace xfdtd {

static std::mutex cout_mutex;

Domain::Domain(std::size_t id, IndexTask task,
               std::shared_ptr<GridSpace> grid_space,
               std::shared_ptr<CalculationParam> calculation_param,
               std::shared_ptr<EMF> emf, std::unique_ptr<Updator> updator,
               std::vector<std::shared_ptr<WaveformSource>> waveform_sources,
               std::vector<std::unique_ptr<Corrector>> correctors,
               std::vector<std::shared_ptr<Monitor>> monitors,
               std::vector<std::shared_ptr<NFFFT>> nfffts,
               std::barrier<>& barrier, bool master)
    : _id{id},
      _task{task},
      _grid_space{std::move(grid_space)},
      _calculation_param{std::move(calculation_param)},
      _emf{std::move(emf)},
      _updator{std::move(updator)},
      _waveform_sources{std::move(waveform_sources)},
      _correctors{std::move(correctors)},
      _monitors{std::move(monitors)},
      _nfffts{std::move(nfffts)},
      _barrier{barrier},
      _master{master} {}

void Domain::run() {
  while (!isCalculationDone()) {
    if (isMaster()) {
      for (auto&& ws : _waveform_sources) {
        ws->updateWaveformSource();
      }
    }

    updateH();

    threadSynchronize();

    correctH();

    synchronize();

    exchangeH();

    synchronize();

    updateE();

    threadSynchronize();

    correctE();

    threadSynchronize();

    record();

    synchronize();

    nextStep();

    threadSynchronize();
  }
}

bool Domain::isCalculationDone() const {
  return _calculation_param->timeParam()->endTimeStep() <=
         _calculation_param->timeParam()->currentTimeStep();
}

void Domain::updateE() { _updator->updateE(); }

void Domain::updateH() { _updator->updateH(); }

void Domain::correctE() {
  for (auto&& c : _correctors) {
    c->correctE();
  }
}

void Domain::correctH() {
  for (auto&& c : _correctors) {
    c->correctH();
  }
}

void Domain::threadSynchronize() { _barrier.arrive_and_wait(); }

void Domain::processSynchronize() {
  if (!isMaster()) {
    return;
  }

  MpiSupport::instance().barrier();
}

void Domain::synchronize() {
  processSynchronize();
  threadSynchronize();
}

void Domain::record() {
  for (auto&& m : _monitors) {
    m->update();
  }

  for (auto&& n : _nfffts) {
    n->update();
  }
}

void Domain::nextStep() {
  if (isMaster() && MpiSupport::instance().isRoot()) {
    std::stringstream ss;
    ss << "\r"
       << "Progress: " << _calculation_param->timeParam()->currentTimeStep() + 1
       << "/" << _calculation_param->timeParam()->endTimeStep();
    std::cout << ss.str() << std::flush;
  }

  if (isMaster()) {
    _calculation_param->timeParam()->nextStep();
  }
}

auto Domain::toString() const -> std::string {
  std::stringstream ss;
  ss << "Domain:\n";
  ss << " ID: " << _id << "\n";
  ss << " Master: " << std::boolalpha << _master << "\n";
  ss << " " << _task.toString() << "\n";
  ss << " " << _updator->toString() << "\n";

  for (auto&& c : _correctors) {
    ss << " " << c->toString() << "\n";
  }

  return ss.str();
}

void Domain::exchangeH() {
  if (!isMaster()) {
    return;
  }

  auto& mpi_support = MpiSupport::instance();

  auto& hx = _emf->hx();
  auto& hy = _emf->hy();
  auto& hz = _emf->hz();

  if (!nodeContainXNBoundary()) {
    mpi_support.sendRecvHyXHead(hy);
    mpi_support.sendRecvHzXHead(hz);
  }

  if (!nodeContainXPBoundary()) {
    mpi_support.recvSendHyXTail(hy);
    mpi_support.recvSendHzXTail(hz);
  }

  if (!nodeContainYNBoundary()) {
    mpi_support.sendRecvHzYHead(hz);
    mpi_support.sendRecvHxYHead(hx);
  }

  if (!nodeContainYPBoundary()) {
    mpi_support.recvSendHzYTail(hz);
    mpi_support.recvSendHxYTail(hx);
  }

  if (!nodeContainZNBoundary()) {
    mpi_support.sendRecvHxZHead(hx);
    mpi_support.sendRecvHyZHead(hy);
  }

  if (!nodeContainZPBoundary()) {
    mpi_support.recvSendHxZTail(hx);
    mpi_support.recvSendHyZTail(hy);
  }

  mpi_support.waitAll();
}

}  // namespace xfdtd
