
#include <memory>

#include "xfdtd/boundary/pml.h"
#include "xfdtd/coordinate_system/coordinate_system.h"
#include "xfdtd/electromagnetic_field/electromagnetic_field.h"
#include "xfdtd/material/material.h"
#include "xfdtd/monitor/current_monitor.h"
#include "xfdtd/monitor/field_monitor.h"
#include "xfdtd/monitor/movie_monitor.h"
#include "xfdtd/monitor/voltage_monitor.h"
#include "xfdtd/network/network.h"
#include "xfdtd/network/port.h"
#include "xfdtd/object/lumped_element/pec_plane.h"
#include "xfdtd/object/lumped_element/resistor.h"
#include "xfdtd/object/lumped_element/voltage_source.h"
#include "xfdtd/object/object.h"
#include "xfdtd/shape/cube.h"
#include "xfdtd/simulation/simulation.h"
#include "xfdtd/waveform/waveform.h"

void microstripLowPass() {
  constexpr double dx{0.4064e-3};
  constexpr double dy{0.4233e-3};
  constexpr double dz{0.265e-3};

  auto domain{std::make_shared<xfdtd::Object>(
      "domain",
      std::make_unique<xfdtd::Cube>(xfdtd::Vector{-5 * dx, -5 * dy, 0 * dz},
                                    xfdtd::Vector{60 * dx, 56 * dy, 13 * dz}),
      xfdtd::Material::createAir())};

  auto substrate{std::make_shared<xfdtd::Object>(
      "substrate",
      std::make_unique<xfdtd::Cube>(xfdtd::Vector{0, 0, 0},
                                    xfdtd::Vector{50 * dx, 46 * dy, 3 * dz}),
      std::make_unique<xfdtd::Material>(
          "sub_material", xfdtd::ElectroMagneticProperty{2.2, 1, 0, 0}))};

  constexpr double center_f{8e9};
  constexpr double band{16e9};
  constexpr double tau{0.996 / band};
  constexpr double t_0{4.5 * tau};

  auto microstrip_0{std::make_shared<xfdtd::PecPlane>(
      "microstrip_0",
      std::make_unique<xfdtd::Cube>(xfdtd::Vector{14 * dx, 0 * dy, 3 * dz},
                                    xfdtd::Vector{6 * dx, 20 * dy, 0 * dz}))};

  auto microstrip_1{std::make_shared<xfdtd::PecPlane>(
      "microstrip_1",
      std::make_unique<xfdtd::Cube>(xfdtd::Vector{30 * dx, 26 * dy, 3 * dz},
                                    xfdtd::Vector{6 * dx, 20 * dy, 0 * dz}))};

  auto microstrip_2{std::make_shared<xfdtd::PecPlane>(
      "microstrip_2",
      std::make_unique<xfdtd::Cube>(xfdtd::Vector{0 * dx, 20 * dy, 3 * dz},
                                    xfdtd::Vector{50 * dx, 6 * dy, 0 * dz}))};

  auto microstrip_3{std::make_shared<xfdtd::PecPlane>(
      "microstrip_3",
      std::make_unique<xfdtd::Cube>(xfdtd::Vector{0 * dx, 0 * dy, 0 * dz},
                                    xfdtd::Vector{50 * dx, 46 * dy, 0 * dz}))};

  auto v_source{std::make_shared<xfdtd::VoltageSource>(
      "voltage_source",
      std::make_unique<xfdtd::Cube>(xfdtd::Vector{14 * dx, 0 * dy, 0 * dz},
                                    xfdtd::Vector{6 * dx, 0 * dy, 3 * dz}),
      xfdtd::Axis::Direction::ZP, 50,
      std::make_unique<xfdtd::Waveform>(
          xfdtd::Waveform::cosineModulatedGaussian(tau, t_0, center_f)))};

  auto r{std::make_shared<xfdtd::Resistor>(
      "resistor",
      std::make_unique<xfdtd::Cube>(xfdtd::Vector{30 * dx, 43 * dy, 0 * dz},
                                    xfdtd::Vector{6 * dx, 0 * dy, 3 * dz}),
      xfdtd::Axis::XYZ::Z, 50)};

  auto c1{std::make_shared<xfdtd::CurrentMonitor>(
      "c1",
      std::make_unique<xfdtd::Cube>(xfdtd::Vector{14 * dx, 10 * dy, 3 * dz},
                                    xfdtd::Vector{6 * dx, 0 * dy, 0 * dz}),
      xfdtd::Axis::Direction::YP, "")};

  auto v1{std::make_shared<xfdtd::VoltageMonitor>(
      "v1",
      std::make_unique<xfdtd::Cube>(xfdtd::Vector{14 * dx, 10 * dy, 0 * dz},
                                    xfdtd::Vector{6 * dx, 0 * dy, 3 * dz}),
      xfdtd::Axis::Direction::ZP, "")};

  auto port_1{std::make_shared<xfdtd::Port>(1, true, 50, c1, v1)};

  auto c2{std::make_shared<xfdtd::CurrentMonitor>(
      "c2",
      std::make_unique<xfdtd::Cube>(xfdtd::Vector{30 * dx, 36 * dy, 3 * dz},
                                    xfdtd::Vector{6 * dx, 0 * dy, 0 * dz}),
      xfdtd::Axis::Direction::YN, "")};

  auto v2{std::make_shared<xfdtd::VoltageMonitor>(
      "v2",
      std::make_unique<xfdtd::Cube>(xfdtd::Vector{30 * dx, 36 * dy, 0 * dz},
                                    xfdtd::Vector{6 * dx, 0 * dy, 3 * dz}),
      xfdtd::Axis::Direction::ZP, "")};

  auto port_2{std::make_shared<xfdtd::Port>(2, false, 50, c2, v2)};

  auto network{std::make_shared<xfdtd::Network>("./data/microstrip_low_pass")};
  network->setFrequencies(xt::arange(2e7, 2e10, 2e7));
  network->addPort(port_1);
  network->addPort(port_2);

  auto ez_monitor{std::make_unique<xfdtd::FieldMonitor>(
      std::make_unique<xfdtd::Cube>(xfdtd::Vector{-5 * dx, -5 * dy, 0},
                                    xfdtd::Vector{60 * dx, 56 * dy, 1 * dz}),
      xfdtd::Axis::XYZ::Z, xfdtd::EMF::Field::EZ)};
  auto movie_monitor{std::make_shared<xfdtd::MovieMonitor>(
      std::move(ez_monitor), 50, "movie", "./data/microstrip_low_pass/movie")};

  auto s{xfdtd::Simulation{dx, dy, dz, 0.90}};
  s.addObject(domain);
  s.addObject(substrate);
  s.addObject(microstrip_0);
  s.addObject(microstrip_1);
  s.addObject(microstrip_2);
  s.addObject(microstrip_3);
  s.addObject(v_source);
  s.addObject(r);
  s.addMonitor(c1);
  s.addMonitor(v1);
  s.addMonitor(c2);
  s.addMonitor(v2);
  s.addNetwork(network);
  s.addBoundary(std::make_shared<xfdtd::PML>(8, xfdtd::Axis::Direction::XN));
  s.addBoundary(std::make_shared<xfdtd::PML>(8, xfdtd::Axis::Direction::XP));
  s.addBoundary(std::make_shared<xfdtd::PML>(8, xfdtd::Axis::Direction::YN));
  s.addBoundary(std::make_shared<xfdtd::PML>(8, xfdtd::Axis::Direction::YP));
  s.addBoundary(std::make_shared<xfdtd::PML>(8, xfdtd::Axis::Direction::ZN));
  s.addBoundary(std::make_shared<xfdtd::PML>(8, xfdtd::Axis::Direction::ZP));
//   s.addMonitor(movie_monitor);

  s.run(3500);

  std::cout << "Output Data..."
            << "\n";
  network->output();
  v2->setOutputDir("./data/microstrip_low_pass/");
  v2->output();
  v1->setOutputDir("./data/microstrip_low_pass/");
  v1->output();
}

int main() {
  // run 20000 time step if don't use PML
  microstripLowPass();
  return 0;
}
