#include <memory>

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
#include "xfdtd/shape/cube.h"
#include "xfdtd/simulation/simulation.h"
#include "xfdtd/waveform/waveform.h"

void microstripBranchLineCoupler() {
  constexpr double dx{0.406e-3};
  constexpr double dy{0.406e-3};
  constexpr double dz{0.265e-3};

  auto domain{std::make_shared<xfdtd::Object>(
      "domain",
      std::make_unique<xfdtd::Cube>(xfdtd::Vector{-25 * dx, -30 * dy, 5 * dz},
                                    xfdtd::Vector{50 * dx, 60 * dy, 13 * dz}),
      xfdtd::Material::createAir())};

  auto substrate{std::make_shared<xfdtd::Object>(
      "substrate",
      std::make_unique<xfdtd::Cube>(xfdtd::Vector{-20 * dx, -25 * dy, 0},
                                    xfdtd::Vector{40 * dx, 50 * dy, 3 * dz}),
      std::make_unique<xfdtd::Material>(
          "sub", xfdtd::ElectroMagneticProperty{2.2, 1, 0, 0}))};

  auto microstrip_1{std::make_shared<xfdtd::PecPlane>(
      "microstrip_1",
      std::make_unique<xfdtd::Cube>(xfdtd::Vector{-12 * dx, -15 * dy, 3 * dz},
                                    xfdtd::Vector{24 * dx, 6 * dy, 0}))};

  auto microstrip_2{std::make_shared<xfdtd::PecPlane>(
      "microstrip_2",
      std::make_unique<xfdtd::Cube>(xfdtd::Vector{-12 * dx, 9 * dy, 3 * dz},
                                    xfdtd::Vector{24 * dx, 6 * dy, 0}))};

  auto microstrip_3{std::make_shared<xfdtd::PecPlane>(
      "microstrip_3",
      std::make_unique<xfdtd::Cube>(xfdtd::Vector{-17 * dx, -12 * dy, 3 * dz},
                                    xfdtd::Vector{10 * dx, 24 * dy, 0}))};

  auto microstrip_4{std::make_shared<xfdtd::PecPlane>(
      "microstrip_4",
      std::make_unique<xfdtd::Cube>(xfdtd::Vector{7 * dx, -12 * dy, 3 * dz},
                                    xfdtd::Vector{10 * dx, 24 * dy, 0}))};

  auto microstrip_5{std::make_shared<xfdtd::PecPlane>(
      "microstrip_5",
      std::make_unique<xfdtd::Cube>(xfdtd::Vector{-15 * dx, -25 * dy, 3 * dz},
                                    xfdtd::Vector{6 * dx, 13 * dy, 0}))};

  auto microstrip_6{std::make_shared<xfdtd::PecPlane>(
      "microstrip_6",
      std::make_unique<xfdtd::Cube>(xfdtd::Vector{9 * dx, -25 * dy, 3 * dz},
                                    xfdtd::Vector{6 * dx, 13 * dy, 0}))};

  auto microstrip_7{std::make_shared<xfdtd::PecPlane>(
      "microstrip_7",
      std::make_unique<xfdtd::Cube>(xfdtd::Vector{-15 * dx, 12 * dy, 3 * dz},
                                    xfdtd::Vector{6 * dx, 13 * dy, 0}))};

  auto microstrip_8{std::make_shared<xfdtd::PecPlane>(
      "microstrip_8",
      std::make_unique<xfdtd::Cube>(xfdtd::Vector{9 * dx, 12 * dy, 3 * dz},
                                    xfdtd::Vector{6 * dx, 13 * dy, 0}))};

  auto ground{std::make_shared<xfdtd::PecPlane>(
      "ground",
      std::make_unique<xfdtd::Cube>(xfdtd::Vector{-20 * dx, -25 * dy, 0},
                                    xfdtd::Vector{40 * dx, 50 * dy, 0}))};

  auto l_min{dz * 30};
  auto tau{l_min / 6e8};
  auto t_0{4.5 * tau};

  auto v_source{std::make_shared<xfdtd::VoltageSource>(
      "voltage_source",
      std::make_unique<xfdtd::Cube>(xfdtd::Vector{-15 * dx, -25 * dy, 0},
                                    xfdtd::Vector{6 * dx, 0 * dy, 3 * dz}),
      xfdtd::Axis::Direction::ZP, 50,
      std::make_unique<xfdtd::Waveform>(xfdtd::Waveform::gaussian(tau, t_0)))};

  auto r1{std::make_shared<xfdtd::Resistor>(
      "resistor_1",
      std::make_unique<xfdtd::Cube>(xfdtd::Vector{9 * dx, -25 * dy, 0},
                                    xfdtd::Vector{6 * dx, 0 * dy, 3 * dz}),
      xfdtd::Axis::XYZ::Z, 50)};

  auto r2{std::make_shared<xfdtd::Resistor>(
      "resistor_2",
      std::make_unique<xfdtd::Cube>(xfdtd::Vector{9 * dx, 25 * dy, 0},
                                    xfdtd::Vector{6 * dx, 0 * dy, 3 * dz}),
      xfdtd::Axis::XYZ::Z, 50)};

  auto r3{std::make_shared<xfdtd::Resistor>(
      "resistor_3",
      std::make_unique<xfdtd::Cube>(xfdtd::Vector{-15 * dx, 25 * dy, 0},
                                    xfdtd::Vector{6 * dx, 0 * dy, 3 * dz}),
      xfdtd::Axis::XYZ::Z, 50)};

  auto v1{std::make_shared<xfdtd::VoltageMonitor>(
      "v1",
      std::make_unique<xfdtd::Cube>(xfdtd::Vector{-15 * dx, -25 * dy, 0},
                                    xfdtd::Vector{6 * dx, 0 * dy, 3 * dz}),
      xfdtd::Axis::Direction::ZP, "")};

  auto v2{std::make_shared<xfdtd::VoltageMonitor>(
      "v2",
      std::make_unique<xfdtd::Cube>(xfdtd::Vector{9 * dx, -25 * dy, 0},
                                    xfdtd::Vector{6 * dx, 0 * dy, 3 * dz}),
      xfdtd::Axis::Direction::ZP, "")};

  auto v3{std::make_shared<xfdtd::VoltageMonitor>(
      "v3",
      std::make_unique<xfdtd::Cube>(xfdtd::Vector{9 * dx, 25 * dy, 0},
                                    xfdtd::Vector{6 * dx, 0 * dy, 3 * dz}),
      xfdtd::Axis::Direction::ZP, "")};

  auto v4{std::make_shared<xfdtd::VoltageMonitor>(
      "v4",
      std::make_unique<xfdtd::Cube>(xfdtd::Vector{-15 * dx, 25 * dy, 0},
                                    xfdtd::Vector{6 * dx, 0 * dy, 3 * dz}),
      xfdtd::Axis::Direction::ZP, "")};

  auto c1{std::make_shared<xfdtd::CurrentMonitor>(
      "c1",
      std::make_unique<xfdtd::Cube>(xfdtd::Vector{-15 * dx, -25 * dy, 2 * dz},
                                    xfdtd::Vector{6 * dx, 0 * dy, 0 * dz}),
      xfdtd::Axis::Direction::ZP, "")};

  auto c2{std::make_shared<xfdtd::CurrentMonitor>(
      "c2",
      std::make_unique<xfdtd::Cube>(xfdtd::Vector{9 * dx, -25 * dy, 2 * dz},
                                    xfdtd::Vector{6 * dx, 0 * dy, 0 * dz}),
      xfdtd::Axis::Direction::ZP, "")};

  auto c3{std::make_shared<xfdtd::CurrentMonitor>(
      "c3",
      std::make_unique<xfdtd::Cube>(xfdtd::Vector{9 * dx, 25 * dy, 2 * dz},
                                    xfdtd::Vector{6 * dx, 0 * dy, 0 * dz}),
      xfdtd::Axis::Direction::ZP, "")};

  auto c4{std::make_shared<xfdtd::CurrentMonitor>(
      "c4",
      std::make_unique<xfdtd::Cube>(xfdtd::Vector{-15 * dx, 25 * dy, 2 * dz},
                                    xfdtd::Vector{6 * dx, 0 * dy, 0 * dz}),
      xfdtd::Axis::Direction::ZP, "")};

  auto port_1{std::make_shared<xfdtd::Port>(1, true, 50, c1, v1)};
  auto port_2{std::make_shared<xfdtd::Port>(2, false, 50, c2, v2)};
  auto port_3{std::make_shared<xfdtd::Port>(3, false, 50, c3, v3)};
  auto port_4{std::make_shared<xfdtd::Port>(4, false, 50, c4, v4)};

  auto network{std::make_shared<xfdtd::Network>(
      std::vector<std::shared_ptr<xfdtd::Port>>{port_1, port_2, port_3, port_4},
      xt::linspace<double>(1e9, 10e9, 100),
      "./data/microstrip_branch_line_coupler")};

  auto s{xfdtd::Simulation{dx, dy, dz, 0.9}};
  s.addObject(domain);
  s.addObject(substrate);
  s.addObject(microstrip_1);
  s.addObject(microstrip_2);
  s.addObject(microstrip_3);
  s.addObject(microstrip_4);
  s.addObject(microstrip_5);
  s.addObject(microstrip_6);
  s.addObject(microstrip_7);
  s.addObject(microstrip_8);
  s.addObject(ground);
  s.addObject(v_source);
  s.addObject(r1);
  s.addObject(r2);
  s.addObject(r3);
  s.addMonitor(v1);
  s.addMonitor(v2);
  s.addMonitor(v3);
  s.addMonitor(v4);
  s.addMonitor(c1);
  s.addMonitor(c2);
  s.addMonitor(c3);
  s.addMonitor(c4);
  s.addNetwork(network);
  s.addBoundary(std::make_shared<xfdtd::PML>(8, xfdtd::Axis::Direction::XN));
  s.addBoundary(std::make_shared<xfdtd::PML>(8, xfdtd::Axis::Direction::XP));
  s.addBoundary(std::make_shared<xfdtd::PML>(8, xfdtd::Axis::Direction::YN));
  s.addBoundary(std::make_shared<xfdtd::PML>(8, xfdtd::Axis::Direction::YP));
  s.addBoundary(std::make_shared<xfdtd::PML>(8, xfdtd::Axis::Direction::ZP));
  s.addBoundary(std::make_shared<xfdtd::PML>(8, xfdtd::Axis::Direction::ZN));
  s.run(4000);

  network->output();
}

int main() {
  microstripBranchLineCoupler();
  return 0;
}