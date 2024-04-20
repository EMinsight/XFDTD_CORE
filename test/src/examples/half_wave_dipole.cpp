
#include <memory>

#include "xfdtd/boundary/pml.h"
#include "xfdtd/common/constant.h"
#include "xfdtd/coordinate_system/coordinate_system.h"
#include "xfdtd/material/material.h"
#include "xfdtd/monitor/current_monitor.h"
#include "xfdtd/monitor/voltage_monitor.h"
#include "xfdtd/network/network.h"
#include "xfdtd/network/port.h"
#include "xfdtd/nffft/nffft.h"
#include "xfdtd/object/lumped_element/voltage_source.h"
#include "xfdtd/object/object.h"
#include "xfdtd/object/thin_wire.h"
#include "xfdtd/shape/cube.h"
#include "xfdtd/shape/cylinder.h"
#include "xfdtd/simulation/simulation.h"

void halfWaveDipole() {
  constexpr double dl{2.5e-4};
  auto domain{std::make_shared<xfdtd::Object>(
      "domain",
      std::make_unique<xfdtd::Cube>(
          xfdtd::Vector(-10 * dl, -10 * dl, -10e-3 - 10 * dl),
          xfdtd::Vector{20 * dl, 20 * dl, 20e-3 + 20 * dl}),
      xfdtd::Material::createAir())};

  auto arm_0{std::make_shared<xfdtd::ThinWire>(
      "arm_0",
      std::make_unique<xfdtd::Cylinder>(xfdtd::Vector{0, 0, 5.125e-3}, 5e-5,
                                        39 * dl, xfdtd::Axis::ZP))};
  auto arm_1{std::make_shared<xfdtd::ThinWire>(
      "arm_1",
      std::make_unique<xfdtd::Cylinder>(xfdtd::Vector{0, 0, -5.125e-3}, 5e-5,
                                        39 * dl, xfdtd::Axis::ZP))};
  constexpr double l_min{dl * 30};
  constexpr double tau{l_min / 6e8};
  constexpr double t_0{4.5 * tau};
  auto v_s{std::make_shared<xfdtd::VoltageSource>(
      "v_s",
      std::make_unique<xfdtd::Cube>(xfdtd::Vector{0, 0, -dl},
                                    xfdtd::Vector{0, 0, 2 * dl}),
      xfdtd::Axis::Direction::ZP, 50, xfdtd::Waveform::gaussian(tau, t_0))};

  auto v1{std::make_shared<xfdtd::VoltageMonitor>(
      "v1",
      std::make_unique<xfdtd::Cube>(xfdtd::Vector{0, 0, -dl},
                                    xfdtd::Vector{0, 0, 2 * dl}),
      xfdtd::Axis::Direction::ZP, "./data/half_wave_dipole")};
  auto c1{std::make_shared<xfdtd::CurrentMonitor>(
      "c1",
      std::make_unique<xfdtd::Cube>(xfdtd::Vector{0, 0, 0},
                                    xfdtd::Vector{0, 0, 0}),
      xfdtd::Axis::Direction::ZP, "./data/half_wave_dipole")};

  auto port{std::make_shared<xfdtd::Port>(1, true, 50, c1, v1)};
  auto network{std::make_shared<xfdtd::Network>(
      std::vector<std::shared_ptr<xfdtd::Port>>{port},
      xt::arange<double>(2e7, 2e10, 2e7), "./data/half_wave_dipole")};

  auto nffft_fd{std::make_shared<xfdtd::NFFFT>(
      13, 13, 13, xt::xarray<double>{7e9}, "./data/half_wave_dipole")};

  auto s{xfdtd::Simulation{dl, dl, dl, 0.9}};
  s.addObject(domain);
  s.addObject(arm_0);
  s.addObject(arm_1);
  s.addObject(v_s);
  s.addMonitor(v1);
  s.addMonitor(c1);
  s.addNetwork(network);
  s.addNF2FF(nffft_fd);
  s.addBoundary(std::make_shared<xfdtd::PML>(8, xfdtd::Axis::Direction::XN));
  s.addBoundary(std::make_shared<xfdtd::PML>(8, xfdtd::Axis::Direction::XP));
  s.addBoundary(std::make_shared<xfdtd::PML>(8, xfdtd::Axis::Direction::YN));
  s.addBoundary(std::make_shared<xfdtd::PML>(8, xfdtd::Axis::Direction::YP));
  s.addBoundary(std::make_shared<xfdtd::PML>(8, xfdtd::Axis::Direction::ZN));
  s.addBoundary(std::make_shared<xfdtd::PML>(8, xfdtd::Axis::Direction::ZP));

  s.run(4000);

  network->output();
  c1->output();
  v1->output();

  nffft_fd->outputRadiationPower();
  nffft_fd->processFarField(
      xfdtd::constant::PI * 0.5,
      xt::linspace<double>(-xfdtd::constant::PI, xfdtd::constant::PI, 360),
      "xy");
  nffft_fd->processFarField(
      xt::linspace<double>(-xfdtd::constant::PI, xfdtd::constant::PI, 360), 0,
      "xz");
  nffft_fd->processFarField(
      xt::linspace<double>(-xfdtd::constant::PI, xfdtd::constant::PI, 360),
      xfdtd::constant::PI * 0.5, "yz");
}

int main() {
  halfWaveDipole();
  return 0;
}
