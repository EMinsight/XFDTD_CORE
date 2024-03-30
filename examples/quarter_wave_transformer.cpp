#include <xfdtd/parallel/parallelized_config.h>

#include <memory>
#include <xtensor.hpp>
#include <xtensor/xnpy.hpp>

#include "xfdtd/boundary/pml.h"
#include "xfdtd/coordinate_system/coordinate_system.h"
#include "xfdtd/material/material.h"
#include "xfdtd/monitor/current_monitor.h"
#include "xfdtd/monitor/voltage_monitor.h"
#include "xfdtd/network/network.h"
#include "xfdtd/network/port.h"
#include "xfdtd/object/lumped_element/pec_plane.h"
#include "xfdtd/object/lumped_element/resistor.h"
#include "xfdtd/object/lumped_element/voltage_source.h"
#include "xfdtd/object/object.h"
#include "xfdtd/parallel/mpi_support.h"
#include "xfdtd/simulation/simulation.h"

void quarterWaveTransformer() {
  constexpr double size{2e-4};

  auto domain{std::make_shared<xfdtd::Object>(
      "domain",
      std::make_unique<xfdtd::Cube>(xfdtd::Vector{0e-3, -1e-3, 0},
                                    xfdtd::Vector{10e-3, 20e-3, 3e-3}),
      xfdtd::Material::createAir())};

  auto substrate{std::make_shared<xfdtd::Object>(
      "substrate",
      std::make_unique<xfdtd::Cube>(xfdtd::Vector{0, 0, 0},
                                    xfdtd::Vector{10e-3, 18e-3, 1e-3}),
      std::make_unique<xfdtd::Material>(
          "sub_material", xfdtd::ElectroMagneticProperty{4.6, 1, 0, 0}))};

  auto microstrip_0{std::make_shared<xfdtd::PecPlane>(
      "microstrip_0",
      std::make_unique<xfdtd::Cube>(xfdtd::Vector{4e-3, 0 * size, 1e-3},
                                    xfdtd::Vector{1.8e-3, 4e-3, 0}))};

  auto microstrip_1{std::make_shared<xfdtd::PecPlane>(
      "microstrip_1",
      std::make_unique<xfdtd::Cube>(xfdtd::Vector{4.4e-3, 4e-3, 1e-3},
                                    xfdtd::Vector{1e-3, 10e-3, 0}))};

  auto microstrip_2{std::make_shared<xfdtd::PecPlane>(
      "microstrip_2",
      std::make_unique<xfdtd::Cube>(xfdtd::Vector{4.8e-3, 14e-3, 1e-3},
                                    xfdtd::Vector{0.4e-3, 4e-3, 0}))};

  auto ground_0{std::make_shared<xfdtd::PecPlane>(
      "ground_0",
      std::make_unique<xfdtd::Cube>(xfdtd::Vector{0 * size, 0, 0},
                                    xfdtd::Vector{10e-3, 18e-3, 0}))};

  constexpr double l_min{size * 30};
  auto tau{l_min / 6e8};
  auto t_0{4.5 * tau};
  auto v_source{std::make_shared<xfdtd::VoltageSource>(
      "voltage_source",
      std::make_unique<xfdtd::Cube>(xfdtd::Vector{4e-3, 0 * size, 0},
                                    xfdtd::Vector{1.8e-3, 0.4e-3, 1e-3}),
      xfdtd::Axis::Direction::ZP, 50,
      std::make_unique<xfdtd::Waveform>(xfdtd::Waveform::gaussian(tau, t_0)))};

  auto r{std::make_shared<xfdtd::Resistor>(
      "resistor",
      std::make_unique<xfdtd::Cube>(xfdtd::Vector{4.8e-3, 17.6e-3, 0},
                                    xfdtd::Vector{0.4e-3, 0.4e-3, 1e-3}),
      xfdtd::Axis::XYZ::Z, 100)};

  auto abc_material{xfdtd::Material(
      "abc_material", xfdtd::ElectroMagneticProperty{1, 142130, 0, 0})};
  auto abs_object_zp{std::make_shared<xfdtd::Object>(
      "abc_object_zp",
      std::make_unique<xfdtd::Cube>(xfdtd::Vector{-1e-3, -2e-3, 3e-3},
                                    xfdtd::Vector{12e-3, 22e-3, 1e-3}),
      std::make_unique<xfdtd::Material>(abc_material))};

  auto abs_object_xn{std::make_shared<xfdtd::Object>(
      "abc_object_xn",
      std::make_unique<xfdtd::Cube>(xfdtd::Vector{-1e-3, -2e-3, 0},
                                    xfdtd::Vector{1e-3, 22e-3, 4e-3}),
      std::make_unique<xfdtd::Material>(abc_material))};

  auto abs_object_xp{std::make_shared<xfdtd::Object>(
      "abc_object_xp",
      std::make_unique<xfdtd::Cube>(xfdtd::Vector{10e-3, -2e-3, 0},
                                    xfdtd::Vector{1e-3, 22e-3, 4e-3}),
      std::make_unique<xfdtd::Material>(abc_material))};

  auto abs_object_yn{std::make_shared<xfdtd::Object>(
      "abc_object_yn",
      std::make_unique<xfdtd::Cube>(xfdtd::Vector{-1e-3, -2e-3, 0},
                                    xfdtd::Vector{12e-3, 1e-3, 4e-3}),
      std::make_unique<xfdtd::Material>(abc_material))};

  auto abs_object_yp{std::make_shared<xfdtd::Object>(
      "abc_object_yp",
      std::make_unique<xfdtd::Cube>(xfdtd::Vector{-1e-3, 19e-3, 0},
                                    xfdtd::Vector{12e-3, 1e-3, 4e-3}),
      std::make_unique<xfdtd::Material>(abc_material))};

  auto v1{std::make_shared<xfdtd::VoltageMonitor>(
      "v1",
      std::make_unique<xfdtd::Cube>(xfdtd::Vector{4e-3, 0, 0},
                                    xfdtd::Vector{1.8e-3, 0.4e-3, 1e-3}),
      xfdtd::Axis::Direction::ZP, "")};

  auto v2{std::make_shared<xfdtd::VoltageMonitor>(
      "v2",
      std::make_unique<xfdtd::Cube>(xfdtd::Vector{4.8e-3, 17.6e-3, 0},
                                    xfdtd::Vector{0.4e-3, 0.4e-3, 1e-3}),
      xfdtd::Axis::Direction::ZP, "")};

  auto c1{std::make_shared<xfdtd::CurrentMonitor>(
      "c1",
      std::make_unique<xfdtd::Cube>(xfdtd::Vector{4e-3, 0, 0.4e-3},
                                    xfdtd::Vector{1.8e-3, 0.4e-3, 0.2e-3}),
      xfdtd::Axis::Direction::ZP, "")};

  auto c2{std::make_shared<xfdtd::CurrentMonitor>(
      "c2",
      std::make_unique<xfdtd::Cube>(xfdtd::Vector{4.8e-3, 17.6e-3, 0.4e-3},
                                    xfdtd::Vector{0.4e-3, 0.4e-3, 0.2e-3}),
      xfdtd::Axis::Direction::ZP, "")};

  auto port_1{std::make_shared<xfdtd::Port>(1, true, 50, c1, v1)};
  auto port_2{std::make_shared<xfdtd::Port>(2, false, 100, c2, v2)};
  auto network{
      std::make_shared<xfdtd::Network>("./data/quarter_wave_transformer")};
  network->addPort(port_1);
  network->addPort(port_2);
  network->setFrequencies(xt::arange(2e7, 8e9, 1e7));

  auto s{
      xfdtd::Simulation{size, size, size, 0.9, xfdtd::ThreadConfig{1, 1, 1}}};
  s.addObject(domain);
  s.addObject(substrate);
  s.addObject(microstrip_0);
  s.addObject(microstrip_1);
  s.addObject(microstrip_2);
  s.addObject(ground_0);
  s.addObject(v_source);
  s.addObject(r);
//   s.addObject(abs_object_zp);
//   s.addObject(abs_object_xn);
//   s.addObject(abs_object_xp);
//   s.addObject(abs_object_yp);
//   s.addObject(abs_object_yn);
  s.addMonitor(v1);
  s.addMonitor(v2);
  s.addMonitor(c1);
  s.addMonitor(c2);
  s.addNetwork(network);
  s.addBoundary(std::make_shared<xfdtd::PML>(-8, xfdtd::Axis::Direction::XN));
  s.addBoundary(std::make_shared<xfdtd::PML>(-8, xfdtd::Axis::Direction::XP));
  s.addBoundary(std::make_shared<xfdtd::PML>(8, xfdtd::Axis::Direction::YN));
  s.addBoundary(std::make_shared<xfdtd::PML>(8, xfdtd::Axis::Direction::YP));
  s.addBoundary(std::make_shared<xfdtd::PML>(8, xfdtd::Axis::Direction::ZP));

  s.run(5000);
  if (xfdtd::MpiSupport::instance().isRoot()) {
    std::cout << "Output Data..."
              << "\n";
  }

  network->output();
  auto output_dir{std::string{"./data/quarter_wave_transformer"}};
  v1->setOutputDir(output_dir);
  v2->setOutputDir(output_dir);
  c1->setOutputDir(output_dir);
  c2->setOutputDir(output_dir);
  v1->output();
  v2->output();
  c1->output();
  c2->output();
  xt::dump_npy(output_dir + "/v_source.npy",
               xt::stack(xt::xtuple(v_source->waveform()->time(),
                                    v_source->waveform()->value())));

  if (xfdtd::MpiSupport::instance().isRoot()) {
    std::cout << "Output Data Finished"
              << "\n";
  }
}

int main() { quarterWaveTransformer(); }
