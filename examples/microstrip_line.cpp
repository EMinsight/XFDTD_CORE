
#include <xfdtd/divider/divider.h>

#include <memory>

#include "xfdtd/boundary/pml.h"
#include "xfdtd/coordinate_system/coordinate_system.h"
#include "xfdtd/material/material.h"
#include "xfdtd/monitor/current_monitor.h"
#include "xfdtd/monitor/voltage_monitor.h"
#include "xfdtd/network/network.h"
#include "xfdtd/network/port.h"
#include "xfdtd/object/lumped_element/pec_plane.h"
#include "xfdtd/object/lumped_element/voltage_source.h"
#include "xfdtd/object/object.h"
#include "xfdtd/parallel/parallelized_config.h"
#include "xfdtd/shape/cube.h"
#include "xfdtd/simulation/simulation.h"
#include "xfdtd/waveform/waveform.h"

void microstripLine() {
  constexpr double dx{0.203e-3};
  constexpr double dy{0.203e-3};
  constexpr double dz{0.1325e-3};

  auto domain{std::make_shared<xfdtd::Object>(
      "domain",
      std::make_unique<xfdtd::Cube>(xfdtd::Vector(0, -4 * dy, 0),
                                    xfdtd::Vector(60 * dx, 64 * dy, 16 * dz)),
      xfdtd::Material::createAir())};

  auto substrate{std::make_shared<xfdtd::Object>(
      "substrate",
      std::make_unique<xfdtd::Cube>(xfdtd::Vector(0, 0, 0),
                                    xfdtd::Vector(60 * dx, 60 * dy, 6 * dz)),
      std::make_unique<xfdtd::Material>(
          "sub", xfdtd::ElectroMagneticProperty{2.2, 1, 0, 0}))};

  auto microstrip_0{std::make_shared<xfdtd::PecPlane>(
      "microstrip_0",
      std::make_unique<xfdtd::Cube>(xfdtd::Vector(24 * dx, 0, 6 * dz),
                                    xfdtd::Vector(12 * dx, 60 * dy, 0)))};
  auto ground_0{std::make_shared<xfdtd::PecPlane>(
      "ground_0",
      std::make_unique<xfdtd::Cube>(xfdtd::Vector(0 * dx, 0, 0 * dz),
                                    xfdtd::Vector(60 * dx, 60 * dy, 0)))};

  auto l_min{30 * dz};
  auto tau{l_min / 6e8};
  auto t_0{tau * 4.5};
  auto v_s{std::make_shared<xfdtd::VoltageSource>(
      "v_s",
      std::make_unique<xfdtd::Cube>(xfdtd::Vector(24 * dx, 0, 0 * dz),
                                    xfdtd::Vector(12 * dx, 0 * dy, 6 * dz)),
      xfdtd::Axis::Direction::ZP, 50,
      std::make_unique<xfdtd::Waveform>(xfdtd::Waveform::gaussian(tau, t_0)))};

  auto v1{std::make_shared<xfdtd::VoltageMonitor>(
      "v1",
      std::make_unique<xfdtd::Cube>(xfdtd::Vector(24 * dx, 40 * dy, 0 * dz),
                                    xfdtd::Vector(12 * dx, 0 * dy, 6 * dz)),
      xfdtd::Axis::Direction::ZP, "")};

  auto c1{std::make_shared<xfdtd::CurrentMonitor>(
      "c1",
      std::make_unique<xfdtd::Cube>(xfdtd::Vector(24 * dx, 40 * dy, 6 * dz),
                                    xfdtd::Vector(12 * dx, 0 * dy, 0 * dz)),
      xfdtd::Axis::Direction::YP, "")};

  auto port_1{std::make_shared<xfdtd::Port>(1, true, 50, c1, v1)};

  auto network{std::make_shared<xfdtd::Network>(
      std::vector<std::shared_ptr<xfdtd::Port>>{port_1},
      xt::linspace(2e7, 10e9, 500), "./data/microstrip_line")};

  auto s{xfdtd::Simulation{dx, dy, dz, 0.9, xfdtd::ThreadConfig{1, 1, 1}}};
  s.addObject(domain);
  s.addObject(substrate);
  s.addObject(microstrip_0);
  s.addObject(ground_0);
  s.addObject(v_s);
  s.addMonitor(v1);
  s.addMonitor(c1);
  s.addNetwork(network);
  s.addBoundary(std::make_shared<xfdtd::PML>(-10, xfdtd::Axis::Direction::XN));
  s.addBoundary(std::make_shared<xfdtd::PML>(-10, xfdtd::Axis::Direction::XP));
  s.addBoundary(std::make_shared<xfdtd::PML>(10, xfdtd::Axis::Direction::YN));
  s.addBoundary(std::make_shared<xfdtd::PML>(-10, xfdtd::Axis::Direction::YP));
  s.addBoundary(std::make_shared<xfdtd::PML>(10, xfdtd::Axis::Direction::ZP));
  s.run(1000);

  v1->setOutputDir("./data/microstrip_line");
  c1->setOutputDir("./data/microstrip_line");
  v1->output();
  c1->output();
  network->output();
}

int main(int argc, char *argv[]) {
  microstripLine();
  return 0;
}
