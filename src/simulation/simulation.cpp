#include <xfdtd/simulation/simulation.h>

#include <chrono>
#include <memory>
#include <thread>
#include <vector>
#include <xtensor/xnpy.hpp>

#include "corrector/corrector.h"
#include "divider/divider.h"
#include "domain/domain.h"
#include "grid_space/grid_space_generator.h"
#include "updator/basic_updator.h"
#include "updator/dispersive_material_updator.h"
#include "updator/updator.h"
#include "xfdtd/calculation_param/calculation_param.h"
#include "xfdtd/grid_space/grid_space.h"
#include "xfdtd/material/dispersive_material.h"
#include "xfdtd/monitor/monitor.h"
#include "xfdtd/nffft/nffft.h"
#include "xfdtd/object/lumped_element/pec_plane.h"
#include "xfdtd/waveform_source/waveform_source.h"

namespace xfdtd {

Simulation::Simulation(double dx, double dy, double dz, double cfl,
                       int num_thread, Divider::Type divider_type)
    : _dx{dx},
      _dy{dy},
      _dz{dz},
      _cfl{cfl},
      _num_thread{num_thread},
      _barrier(num_thread),
      _divider_type{divider_type} {}

Simulation::~Simulation() = default;

void Simulation::addObject(std::shared_ptr<xfdtd::Object> object) {
  _objects.emplace_back(std::move(object));
}

void Simulation::addWaveformSource(
    std::shared_ptr<WaveformSource> waveform_source) {
  _waveform_sources.emplace_back(std::move(waveform_source));
}

void Simulation::addBoundary(std::shared_ptr<Boundary> boundary) {
  _boundaries.emplace_back(std::move(boundary));
}

void Simulation::addMonitor(std::shared_ptr<Monitor> monitor) {
  _monitors.emplace_back(std::move(monitor));
}

void Simulation::addNetwork(std::shared_ptr<Network> network) {
  _networks.emplace_back(std::move(network));
}

void Simulation::addNF2FF(std::shared_ptr<NFFFT> nffft) {
  _nfffts.emplace_back(std::move(nffft));
}

const std::shared_ptr<CalculationParam>& Simulation::calculationParam() const {
  return _calculation_param;
}

const std::shared_ptr<GridSpace>& Simulation::gridSpace() const {
  return _grid_space;
}

const std::shared_ptr<EMF>& Simulation::emf() const { return _emf; }

void Simulation::run(std::size_t time_step) {
  _start_time = std::chrono::high_resolution_clock::now();
  std::cout << "Simulation start..."
            << "\n";
  init(time_step);

  struct ThreadGuard {
    ~ThreadGuard() {
      for (auto& t : _threads) {
        if (t.joinable()) {
          t.join();
        }
      }
    }

    std::vector<std::thread> _threads;
    std::chrono::high_resolution_clock::time_point _start_time;
  };

  {
    std::vector<std::thread> threads;
    for (std::size_t i = 1; i < _domains.size(); ++i) {
      threads.emplace_back([&domain = _domains[i]]() { domain->run(); });
    }

    ThreadGuard tg{
        ._threads = std::move(threads),
    };

    _domains[0]->run();
  }

  _end_time = std::chrono::high_resolution_clock::now();
  std::cout << "\n"
            << "Elapsed time: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(
                   _end_time - _start_time)
                   .count()
            << " ms"
            << " or "
            << std::chrono::duration_cast<std::chrono::seconds>(_end_time -
                                                                _start_time)
                   .count()
            << " s"
            << "\n";
}

void Simulation::init(std::size_t time_step) {
  // First: generate grid space
  generateGridSpace();

  // Second: allocate the basic concept: space,param,emf
  _calculation_param =
      std::make_shared<CalculationParam>(_grid_space.get(), 0, time_step, _cfl);
  _emf = std::make_shared<EMF>();

  // Third: init all the objects
  for (const auto& o : _objects) {
    o->init(_grid_space, _calculation_param, _emf);
  }
  for (const auto& b : _boundaries) {
    b->init(_grid_space, _calculation_param, _emf);
  }
  for (const auto& s : _waveform_sources) {
    s->init(_grid_space, _calculation_param, _emf);
  }

  // correct material space: param
  correctMaterialSpace();
  _calculation_param->calculateCoefficient(_grid_space.get());
  correctUpdateCoefficient();

  // init monitor
  for (const auto& m : _monitors) {
    m->init(_grid_space, _calculation_param, _emf);
  }
  for (const auto& n : _networks) {
    n->init(_grid_space, _calculation_param, _emf);
  }
  for (const auto& n : _nfffts) {
    n->init(_grid_space, _calculation_param, _emf);
  }

  _emf->allocateEx(_grid_space->sizeX(), _grid_space->sizeY() + 1,
                   _grid_space->sizeZ() + 1);
  _emf->allocateEy(_grid_space->sizeX() + 1, _grid_space->sizeY(),
                   _grid_space->sizeZ() + 1);
  _emf->allocateEz(_grid_space->sizeX() + 1, _grid_space->sizeY() + 1,
                   _grid_space->sizeZ());
  _emf->allocateHx(_grid_space->sizeX() + 1, _grid_space->sizeY(),
                   _grid_space->sizeZ());
  _emf->allocateHy(_grid_space->sizeX(), _grid_space->sizeY() + 1,
                   _grid_space->sizeZ());
  _emf->allocateHz(_grid_space->sizeX(), _grid_space->sizeY(),
                   _grid_space->sizeZ() + 1);
  generateDomain();
}

void Simulation::generateDomain() {
  if (_num_thread <= 1) {
    _num_thread = 1;
  }

  if (std::thread::hardware_concurrency() < _num_thread) {
    std::cout << "The number of threads is too large, set to the maximum "
                 "number of threads: "
              << std::thread::hardware_concurrency() << "\n";
    _num_thread = std::thread::hardware_concurrency();
  }

  if (_num_thread == 1) {
    std::cout << "Single thread mode\n";
  }

  Divider::IndexTask problem = Divider::makeTask(
      Divider::makeRange<std::size_t>(0, _grid_space->sizeX()),
      Divider::makeRange<std::size_t>(0, _grid_space->sizeY()),
      Divider::makeRange<std::size_t>(0, _grid_space->sizeZ()));

  auto tasks = std::vector<Divider::IndexTask>{Divider::makeTask(
      Divider::makeRange<std::size_t>(0, _grid_space->sizeX()),
      Divider::makeRange<std::size_t>(0, _grid_space->sizeY()),
      Divider::makeRange<std::size_t>(0, _grid_space->sizeZ()))};

  try {
    tasks = Divider::divide(problem, _num_thread, _divider_type);
  } catch (const XFDTDDividerException& e) {
    std::cout << "Single thread mode\n";
  }

  // Corrector

  std::vector<std::unique_ptr<Corrector>> correctors;

  for (auto&& w : _waveform_sources) {
    auto c = w->generateCorrector(problem);
    if (c == nullptr) {
      continue;
    }

    correctors.emplace_back(std::move(c));
  }

  for (auto&& o : _objects) {
    auto c = o->generateCorrector(problem);
    if (c == nullptr) {
      continue;
    }

    correctors.emplace_back(std::move(c));
  }

  for (auto&& b : _boundaries) {
    auto c = b->generateDomainCorrector(problem);
    if (c == nullptr) {
      continue;
    }

    correctors.emplace_back(std::move(c));
  }

  bool master = true;
  std::size_t id = {0};
  for (const auto& t : tasks) {
    auto updator = makeUpdator(t);
    if (master) {
      _domains.emplace_back(std::make_unique<Domain>(
          id, t, _grid_space, _calculation_param, _emf, std::move(updator),
          _waveform_sources, std::move(correctors), _monitors, _nfffts,
          _barrier, master));
      master = false;
    } else {
      _domains.emplace_back(std::make_unique<Domain>(
          id, t, _grid_space, _calculation_param, _emf, std::move(updator),
          std::vector<std::shared_ptr<WaveformSource>>{},
          std::vector<std::unique_ptr<Corrector>>{},
          std::vector<std::shared_ptr<Monitor>>{},
          std::vector<std::shared_ptr<NFFFT>>{}, _barrier, master));
    }
    ++id;
  }
}

void Simulation::generateGridSpace() {
  std::vector<std::shared_ptr<Shape>> shapes;
  shapes.reserve(_objects.size());
  for (const auto& o : _objects) {
    shapes.emplace_back(o->shape()->clone());
  }

  if (!_boundaries.empty()) {
    _grid_space =
        GridSpaceGenerator::generate(shapes, _boundaries, _dx, _dy, _dz);
    return;
  }

  _grid_space = GridSpaceGenerator::generate(shapes, _dx, _dy, _dz);
}

void Simulation::correctMaterialSpace() {
  std::size_t m_index = {0};
  for (auto&& o : _objects) {
    if (std::dynamic_pointer_cast<PecPlane>(o) != nullptr) {
      continue;
    }

    o->correctMaterialSpace(m_index);
    _calculation_param->materialParam()->addMaterial(o->material());
    ++m_index;
  }

  for (auto&& b : _boundaries) {
    b->correctMaterialSpace();
  }

  for (auto&& s : _waveform_sources) {
    s->correctMaterialSpace();
  }

  _calculation_param->generateMaterialSpaceParam(_grid_space.get());
  for (auto&& o : _objects) {
    if (std::dynamic_pointer_cast<PecPlane>(o) == nullptr) {
      continue;
    }
    // for pec plane
    o->correctMaterialSpace(-1);
  }
}

void Simulation::correctUpdateCoefficient() {
  for (auto&& o : _objects) {
    o->correctUpdateCoefficient();
  }

  for (auto&& b : _boundaries) {
    b->correctUpdateCoefficient();
  }

  for (auto&& s : _waveform_sources) {
    s->correctUpdateCoefficient();
  }
}

std::unique_ptr<Updator> Simulation::makeUpdator(
    const Divider::IndexTask& task) {
  bool dispersion = false;
  for (const auto& m : _calculation_param->materialParam()->materialArray()) {
    if (m->dispersion()) {
      dispersion = true;
      break;
    }
  }

  if (!dispersion) {
    if (_grid_space->dimension() == GridSpace::Dimension::THREE) {
      return std::make_unique<BasicUpdator3D>(_grid_space, _calculation_param,
                                              _emf, task);
    }
    if (_grid_space->dimension() == GridSpace::Dimension::TWO) {
      return std::make_unique<BasicUpdatorTE>(_grid_space, _calculation_param,
                                              _emf, task);
    }
    throw XFDTDSimulationException("Invalid dimension");
  }

  // Contains linear dispersive material
  if (dispersion && _grid_space->dimension() == GridSpace::Dimension::THREE) {
    for (const auto& m : _calculation_param->materialParam()->materialArray()) {
      if (!m->dispersion()) {
        continue;
      }

      auto dispersion_material =
          std::dynamic_pointer_cast<LinearDispersiveMaterial>(m);
      if (!dispersion_material) {
        break;
      }

      if (auto lorentz_material =
              std::dynamic_pointer_cast<LorentzMedium>(dispersion_material);
          lorentz_material != nullptr) {
        _emf->allocateExPrev(_grid_space->sizeX(), _grid_space->sizeY() + 1,
                             _grid_space->sizeZ() + 1);
        _emf->allocateEyPrev(_grid_space->sizeX() + 1, _grid_space->sizeY(),
                             _grid_space->sizeZ() + 1);
        _emf->allocateEzPrev(_grid_space->sizeX() + 1, _grid_space->sizeY() + 1,
                             _grid_space->sizeZ());

        return std::make_unique<LorentzADEUpdator>(
            _grid_space, _calculation_param, _emf, task);
      }

      if (auto drude_material =
              std::dynamic_pointer_cast<DrudeMedium>(dispersion_material);
          drude_material != nullptr) {
        return std::make_unique<DrudeADEUpdator>(
            _grid_space, _calculation_param, _emf, task);
      }

      if (auto debye_material =
              std::dynamic_pointer_cast<DebyeMedium>(dispersion_material);
          debye_material != nullptr) {
        return std::make_unique<DebyeADEUpdator>(
            _grid_space, _calculation_param, _emf, task);
      }
    }
  }

  throw XFDTDSimulationException("don't support this type of updator yet.");
}

}  // namespace xfdtd
