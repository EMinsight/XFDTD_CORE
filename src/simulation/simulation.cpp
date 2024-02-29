#include <xfdtd/simulation/simulation.h>

#include <chrono>
#include <memory>
#include <xtensor/xnpy.hpp>

#include "divider/divider.h"
#include "grid_space/grid_space_generator.h"
#include "updator/basic_updator.h"
#include "updator/dispersive_material_updator.h"
#include "updator/updator.h"
#include "xfdtd/calculation_param/calculation_param.h"
#include "xfdtd/grid_space/grid_space.h"
#include "xfdtd/material/dispersive_material.h"
#include "xfdtd/object/lumped_element/pec_plane.h"

namespace xfdtd {

Simulation::Simulation(double dx, double dy, double dz, double cfl)
    : _dx{dx}, _dy{dy}, _dz{dz}, _cfl{cfl}, _time_step{0} {}

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

  // xt::dump_npy("./data/check/cexe.npy",
  //              _calculation_param->fdtdCoefficient()->cexe());
  // xt::dump_npy("./data/check/cexhy.npy",
  //              _calculation_param->fdtdCoefficient()->cexhy());
  // xt::dump_npy("./data/check/cexhz.npy",
  //              _calculation_param->fdtdCoefficient()->cexhz());
  // xt::dump_npy("./data/check/ceye.npy",
  //              _calculation_param->fdtdCoefficient()->ceye());
  // xt::dump_npy("./data/check/ceyhz.npy",
  //              _calculation_param->fdtdCoefficient()->ceyhz());
  // xt::dump_npy("./data/check/ceyhx.npy",
  //              _calculation_param->fdtdCoefficient()->ceyhx());
  // xt::dump_npy("./data/check/ceze.npy",
  //              _calculation_param->fdtdCoefficient()->ceze());
  // xt::dump_npy("./data/check/cezhx.npy",
  //              _calculation_param->fdtdCoefficient()->cezhx());
  // xt::dump_npy("./data/check/cezhy.npy",
  //              _calculation_param->fdtdCoefficient()->cezhy());
  // xt::dump_npy("./data/check/chxh.npy",
  //              _calculation_param->fdtdCoefficient()->chxh());
  // xt::dump_npy("./data/check/chxey.npy",
  //              _calculation_param->fdtdCoefficient()->chxey());
  // xt::dump_npy("./data/check/chxez.npy",
  //              _calculation_param->fdtdCoefficient()->chxez());
  // xt::dump_npy("./data/check/chyh.npy",
  //              _calculation_param->fdtdCoefficient()->chyh());
  // xt::dump_npy("./data/check/chyex.npy",
  //              _calculation_param->fdtdCoefficient()->chyex());
  // xt::dump_npy("./data/check/chyez.npy",
  //              _calculation_param->fdtdCoefficient()->chyez());
  // xt::dump_npy("./data/check/chzh.npy",
  //              _calculation_param->fdtdCoefficient()->chzh());
  // xt::dump_npy("./data/check/chzex.npy",
  //              _calculation_param->fdtdCoefficient()->chzex());
  // xt::dump_npy("./data/check/chzey.npy",
  //              _calculation_param->fdtdCoefficient()->chzey());
  // xt::dump_npy("./data/check/eps_x.npy",
  //              _calculation_param->materialParam()->epsX());
  // xt::dump_npy("./data/check/sigma_e_x.npy",
  //              _calculation_param->materialParam()->sigmaEX());

  // const std::shared_ptr<const GridSpace> grid_space = _grid_space;

  // const auto& h_size_x{grid_space->hSizeX()};
  // const auto& h_size_y{grid_space->hSizeY()};
  // const auto& h_size_z{grid_space->hSizeZ()};
  // const auto& e_size_x{grid_space->eSizeX()};
  // const auto& e_size_y{grid_space->eSizeY()};
  // const auto& e_size_z{grid_space->eSizeZ()};

  // auto e_x_size{xt::meshgrid(e_size_x, h_size_y, h_size_z)};
  // xt::dump_npy("./data/check/dz.npy", std::get<2>(e_x_size));
  // std::cout << _calculation_param->timeParam()->dt() << "\n";
  // exit(0);

  while (_calculation_param->timeParam()->currentTimeStep() <
         _calculation_param->timeParam()->endTimeStep()) {
    updateE();
    correctE();
    updateH();
    correctH();

    record();
    printRunInfo();
    _calculation_param->timeParam()->nextStep();
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
  _time_step = time_step;
  // First: generate grid space
  generateGridSpace();

  // Second: allocate the basic concept: space,param,emf
  _calculation_param = std::make_shared<CalculationParam>(_grid_space.get(), 0,
                                                          _time_step, _cfl);
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

  setUpdator();

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

  for (auto&& o : _objects) {
    // if (auto v = std::dynamic_pointer_cast<Inductor>(o); v != nullptr) {
    //   continue;
    // }

    auto c = o->generateCorrector(
        Divider::Task<std::size_t>{{0, _grid_space->sizeX()},
                                   {0, _grid_space->sizeY()},
                                   {0, _grid_space->sizeZ()}});
    if (c == nullptr) {
      continue;
    }

    _correctors.emplace_back(std::move(c));
  }

  for (auto&& b : _boundaries) {
    auto c = b->generateDomainCorrector(
        Divider::Task<std::size_t>{{0, _grid_space->sizeX()},
                                   {0, _grid_space->sizeY()},
                                   {0, _grid_space->sizeZ()}});
    if (c == nullptr) {
      continue;
    }

    _correctors.emplace_back(std::move(c));
  }
}

void Simulation::updateE() {
  for (const auto& s : _waveform_sources) {
    s->updateWaveformSourceE();
  }
  _updator->updateE();
}

void Simulation::updateH() {
  for (const auto& s : _waveform_sources) {
    s->updateWaveformSourceH();
  }
  _updator->updateH();
}

void Simulation::correctE() {
  for (const auto& s : _waveform_sources) {
    s->correctE();
  }

  for (auto&& c : _correctors) {
    c->correctE();
  }
}

void Simulation::correctH() {
  for (const auto& s : _waveform_sources) {
    s->correctH();
  }

  for (auto&& c : _correctors) {
    c->correctH();
  }
}

void Simulation::record() {
  for (const auto& m : _monitors) {
    m->update();
  }

  for (const auto& n : _nfffts) {
    n->update();
  }
}

void Simulation::printRunInfo() {
  std::cout << "\r";
  std::cout << "Progress: "
            << _calculation_param->timeParam()->currentTimeStep() + 1 << "/"
            << _calculation_param->timeParam()->endTimeStep() << std::flush;
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

void Simulation::setUpdator() {
  bool dispersion = false;
  for (const auto& m : _calculation_param->materialParam()->materialArray()) {
    if (m->dispersion()) {
      dispersion = true;
      break;
    }
  }

  if (!dispersion) {
    if (_grid_space->dimension() == GridSpace::Dimension::THREE) {
      _updator = std::make_unique<BasicUpdator3D>(_grid_space,
                                                  _calculation_param, _emf);
    } else if (_grid_space->dimension() == GridSpace::Dimension::TWO) {
      _updator = std::make_unique<BasicUpdatorTE>(_grid_space,
                                                  _calculation_param, _emf);
    } else {
      throw std::runtime_error("Invalid dimension");
    }
    return;
  }

  // Contains linear dispersive material
  if (dispersion && _grid_space->dimension() == GridSpace::Dimension::THREE) {
    _updator = std::make_unique<LorentzADEUpdator>(_grid_space,
                                                   _calculation_param, _emf);
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
        _updator = std::make_unique<LorentzADEUpdator>(
            _grid_space, _calculation_param, _emf);
        std::cout << "\nDecide to use LorentzADEUpdator\n";
        return;
      }

      if (auto drude_material =
              std::dynamic_pointer_cast<DrudeMedium>(dispersion_material);
          drude_material != nullptr) {
        _updator = std::make_unique<DrudeADEUpdator>(_grid_space,
                                                     _calculation_param, _emf);
        std::cout << "\nDecide to use DrudeADEUpdator\n";
        return;
      }

      if (auto debye_material =
              std::dynamic_pointer_cast<DebyeMedium>(dispersion_material);
          debye_material != nullptr) {
        _updator = std::make_unique<DebyeADEUpdator>(_grid_space,
                                                     _calculation_param, _emf);
        std::cout << "\nDecide to use DebyeADEUpdator\n";
        return;
      }
    }
  }

  throw XFDTDSimulationException("don't support this type of updator yet.");
}

}  // namespace xfdtd
