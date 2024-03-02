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

void Domain::communicateForH() {
  // only support for shared memory model
  _updator->setHyBufferForEdgeZN(makeHyBufferForEx());
  _updator->setHzBufferForEdgeYN(makeHzBufferForEx());

  _updator->setHzBufferForEdgeXN(makeHzBufferForEy());
  _updator->setHxBufferForEdgeZN(makeHxBufferForEy());

  _updator->setHxBufferForEdgeYN(makeHxBufferForEz());
  _updator->setHyBufferForEdgeXN(makeHyBufferForEz());
  synchronize();
}

void Domain::communicateForE() {}

// TODO(franzero): fixed all make buffer functions

xt::xarray<double> Domain::makeHyBufferForEx() {
  if (containZNEdge()) {
    return {};
  }

  const auto x_range = _task.xRange();
  const auto y_range = _task.yRange();
  const auto k = _task.zRange().start();
  xt::xarray<double> hy_buffer =
      xt::zeros<double>({x_range.size(), y_range.size(), std::size_t{1}});
  for (std::size_t i{x_range.start()}; i < x_range.end(); ++i) {
    for (std::size_t j{y_range.start()}; j < y_range.end(); ++j) {
      hy_buffer(i, j, 0) = _emf->hy()(i, j, k - 1);
    }
  }

  return hy_buffer;
}

xt::xarray<double> Domain::makeHzBufferForEx() {
  if (containYNEdge()) {
    return {};
  }

  const auto x_range = _task.xRange();
  const auto z_range = _task.zRange();
  const auto j = _task.yRange().start();
  xt::xarray<double> hz_buffer =
      xt::zeros<double>({x_range.size(), std::size_t{1}, z_range.size()});
  for (std::size_t i{x_range.start()}; i < x_range.end(); ++i) {
    for (std::size_t k{z_range.start()}; k < z_range.end(); ++k) {
      hz_buffer(i, 0, k) = _emf->hz()(i, j - 1, k);
    }
  }

  return hz_buffer;
}

xt::xarray<double> Domain::makeHzBufferForEy() {
  if (containXNEdge()) {
    return {};
  }

  const auto y_range = _task.yRange();
  const auto z_range = _task.zRange();
  const auto i = _task.xRange().start();
  xt::xarray<double> hz_buffer =
      xt::zeros<double>({std::size_t{1}, y_range.size(), z_range.size()});
  for (std::size_t j{y_range.start()}; j < y_range.end(); ++j) {
    for (std::size_t k{z_range.start()}; k < z_range.end(); ++k) {
      hz_buffer(0, j, k) = _emf->hz()(i - 1, j, k);
    }
  }

  return hz_buffer;
}

xt::xarray<double> Domain::makeHxBufferForEy() {
  if (containZNEdge()) {
    return {};
  }

  const auto x_range = _task.xRange();
  const auto y_range = _task.yRange();
  const auto k = _task.zRange().start();
  xt::xarray<double> hx_buffer =
      xt::zeros<double>({x_range.size(), y_range.size(), std::size_t{1}});
  for (std::size_t i{x_range.start()}; i < x_range.end(); ++i) {
    for (std::size_t j{y_range.start()}; j < y_range.end(); ++j) {
      hx_buffer(i, j, 0) = _emf->hx()(i, j, k - 1);
    }
  }

  return hx_buffer;
}

xt::xarray<double> Domain::makeHxBufferForEz() {
  // for Ez
  if (containYNEdge()) {
    return {};
  }

  const auto x_range = _task.xRange();
  const auto z_range = _task.zRange();
  const auto j = _task.yRange().start();
  xt::xarray<double> hx_buffer =
      xt::zeros<double>({x_range.size(), std::size_t{1}, z_range.size()});
  const auto is = x_range.start();
  const auto ie = x_range.end();
  const auto ks = z_range.start();
  const auto ke = z_range.end();
  for (std::size_t i{is}; i < ie; ++i) {
    for (std::size_t k{ks}; k < ke; ++k) {
      hx_buffer(i - is, 0, k - ks) = _emf->hx()(i, j - 1, k);
    }
  }

  return hx_buffer;
}

xt::xarray<double> Domain::makeHyBufferForEz() {
  if (containXNEdge()) {
    return {};
  }

  const auto y_range = _task.yRange();
  const auto z_range = _task.zRange();
  const auto i = _task.xRange().start();
  xt::xarray<double> hy_buffer =
      xt::zeros<double>({std::size_t{1}, y_range.size(), z_range.size()});
  for (std::size_t j{y_range.start()}; j < y_range.end(); ++j) {
    for (std::size_t k{z_range.start()}; k < z_range.end(); ++k) {
      hy_buffer(0, j, k) = _emf->hy()(i - 1, j, k);
    }
  }

  return hy_buffer;
}

}  // namespace xfdtd
