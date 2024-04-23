#include <filesystem>
#include <memory>
#include <xtensor.hpp>
#include <xtensor/xnpy.hpp>

#include "xfdtd/boundary/pml.h"
#include "xfdtd/common/constant.h"
#include "xfdtd/coordinate_system/coordinate_system.h"
#include "xfdtd/material/dispersive_material.h"
#include "xfdtd/material/material.h"
#include "xfdtd/monitor/field_monitor.h"
#include "xfdtd/monitor/field_time_monitor.h"
#include "xfdtd/monitor/movie_monitor.h"
#include "xfdtd/object/object.h"
#include "xfdtd/parallel/parallelized_config.h"
#include "xfdtd/shape/cube.h"
#include "xfdtd/simulation/simulation.h"
#include "xfdtd/waveform/waveform.h"
#include "xfdtd/waveform_source/tfsf_1d.h"

auto slabRefectionTransmission() {
  const std::string output_dir = "./data/slab_reflection_transmission";
  const auto path = std::filesystem::path(output_dir);
  constexpr double dl = 0.5e-3;

  auto domain = std::make_shared<xfdtd::Object>(
      "domain",
      std::make_unique<xfdtd::Cube>(
          xfdtd::Vector{xfdtd::constant::NEG_INF, xfdtd::constant::NEG_INF, 0},
          xfdtd::Vector{xfdtd::constant::INF, xfdtd::constant::INF, 300 * dl}),
      xfdtd::Material::createAir());

  auto slab = std::make_shared<xfdtd::Object>(
      "slab",
      std::make_unique<xfdtd::Cube>(
          xfdtd::Vector{xfdtd::constant::NEG_INF, xfdtd::constant::NEG_INF,
                        50 * dl},
          xfdtd::Vector{xfdtd::constant::INF, xfdtd::constant::INF, 260 * dl}),
      xfdtd::MLorentzMaterial::makeMLorentz("human_blood", 31.1662, {6.9379e21},
                                            {1.5057e12}, {6.1637e18},
                                            {4.5425e10}, {1}));

  auto l_min{dl * 40};
  auto freq_max = 3e8 / (dl * 40);
  auto tau{l_min / 6e8};
  auto t_0{4.5 * tau};

  auto tfsf_1d = std::make_shared<xfdtd::TFSF1D>(
      25, true, xfdtd::Waveform::gaussian(tau, t_0));

  auto reflect_ex = std::make_shared<xfdtd::FieldTimeMonitor>(
      "reflect_ex",
      std::make_unique<xfdtd::Cube>(
          xfdtd::Vector{xfdtd::constant::NEG_INF, xfdtd::constant::NEG_INF,
                        10 * dl},
          xfdtd::Vector{xfdtd::constant::INF, xfdtd::constant::INF, dl}),
      output_dir, xfdtd::EMF::Field::EX);

  auto transmit_ex = std::make_shared<xfdtd::FieldTimeMonitor>(
      "transmit_ex",
      std::make_unique<xfdtd::Cube>(
          xfdtd::Vector{xfdtd::constant::NEG_INF, xfdtd::constant::NEG_INF,
                        50 * dl},
          xfdtd::Vector{xfdtd::constant::INF, xfdtd::constant::INF, dl}),
      output_dir, xfdtd::EMF::Field::EX);

  auto ex_movie = std::make_shared<xfdtd::MovieMonitor>(
      std::make_unique<xfdtd::FieldMonitor>(
          std::make_unique<xfdtd::Cube>(
              xfdtd::Vector{xfdtd::constant::NEG_INF, xfdtd::constant::NEG_INF,
                            0},
              xfdtd::Vector{xfdtd::constant::INF, xfdtd::constant::INF,
                            200 * dl}),
          xfdtd::EMF::Field::EX, "", ""),
      5, "movie", (path / "ex_movie").string());

  auto s{xfdtd::Simulation{dl, dl, dl, 0.9, xfdtd::ThreadConfig{1, 1, 2}}};
  s.addObject(domain);
  s.addObject(slab);
  s.addWaveformSource(tfsf_1d);
  s.addMonitor(reflect_ex);
  s.addMonitor(transmit_ex);
  s.addMonitor(ex_movie);
  s.addBoundary(std::make_shared<xfdtd::PML>(10, xfdtd::Axis::Direction::ZP));
  s.addBoundary(std::make_shared<xfdtd::PML>(10, xfdtd::Axis::Direction::ZN));
  s.run(700);

  reflect_ex->output();
  transmit_ex->output();

  if (s.isRoot()) {
    auto time = tfsf_1d->waveform()->time();
    auto value = tfsf_1d->waveform()->value();
    xt::dump_npy((path / "incident.npy").string(),
                 xt::stack(xt::xtuple(time, value)));

    auto human_blood =
        std::dynamic_pointer_cast<xfdtd::MLorentzMaterial>(slab->material());
    if (human_blood == nullptr) {
      throw std::runtime_error("slab material is not MLorentzMaterial");
    }

    auto epsilon = human_blood->relativePermittivity(
        xt::linspace(freq_max / 100, freq_max, 100));
    xt::dump_npy((path / "epsilon.npy").string(),
                 xt::stack(xt::xtuple(
                     xt::linspace(freq_max / 100, freq_max, 100), epsilon)));
  }
}

int main(int argc, char const *argv[]) {
  slabRefectionTransmission();
  return 0;
}
