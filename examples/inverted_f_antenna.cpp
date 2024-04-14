#include <memory>

#include "xfdtd/boundary/pml.h"
#include "xfdtd/common/constant.h"
#include "xfdtd/coordinate_system/coordinate_system.h"
#include "xfdtd/material/material.h"
#include "xfdtd/monitor/current_monitor.h"
#include "xfdtd/monitor/voltage_monitor.h"
#include "xfdtd/network/network.h"
#include "xfdtd/network/port.h"
#include "xfdtd/object/lumped_element/pec_plane.h"
#include "xfdtd/object/lumped_element/voltage_source.h"
#include "xfdtd/object/object.h"
#include "xfdtd/shape/cube.h"
#include "xfdtd/simulation/simulation.h"
#include "xfdtd/waveform/waveform.h"

void invertedFAntenna() {
  constexpr double dx{0.262e-3};
  constexpr double dy{0.4e-3};
  constexpr double dz{0.4e-3};

  auto domain{std::make_shared<xfdtd::Object>(
      "domain",
      std::make_unique<xfdtd::Cube>(
          xfdtd::Vector{-0.787e-3 - 10 * dx, -10 * dy, -10 * dz},
          xfdtd::Vector{0.787e-3 + 20 * dx, 40e-3 + 20 * dy, 40e-3 + 20 * dz}),
      xfdtd::Material::createAir())};

  auto substrate{std::make_shared<xfdtd::Object>(
      "substrate",
      std::make_unique<xfdtd::Cube>(xfdtd::Vector{-0.787e-3, 0, 0},
                                    xfdtd::Vector{0.787e-3, 40e-3, 40e-3}),
      std::make_unique<xfdtd::Material>(
          "sub", xfdtd::ElectroMagneticProperty{2.2, 1, 0, 0}))};

  auto antenna_0{std::make_shared<xfdtd::PecPlane>(
      "0", std::make_unique<xfdtd::Cube>(xfdtd::Vector{0, 0, 24e-3},
                                         xfdtd::Vector{0, 28.4e-3, 2.4e-3}))};

  auto antenna_1{std::make_shared<xfdtd::PecPlane>(
      "1", std::make_unique<xfdtd::Cube>(xfdtd::Vector{0, 16e-3, 30e-3},
                                         xfdtd::Vector{0, 12.4e-3, 2.4e-3}))};

  auto antenna_2{std::make_shared<xfdtd::PecPlane>(
      "2", std::make_unique<xfdtd::Cube>(xfdtd::Vector{0, 26e-3, 8.4e-3},
                                         xfdtd::Vector{0, 2.4e-3, 24e-3}))};

  auto antenna_3{std::make_shared<xfdtd::PecPlane>(
      "3", std::make_unique<xfdtd::Cube>(xfdtd::Vector{0, 20.8e-3, 16e-3},
                                         xfdtd::Vector{0, 2.4e-3, 16.4e-3}))};

  auto antenna_4{std::make_shared<xfdtd::PecPlane>(
      "4", std::make_unique<xfdtd::Cube>(xfdtd::Vector{-0.787e-3, 16e-3, 30e-3},
                                         xfdtd::Vector{0.787e-3, 0, 2.4e-3}))};

  auto ground{std::make_shared<xfdtd::PecPlane>(
      "ground", std::make_unique<xfdtd::Cube>(xfdtd::Vector{-0.787e-3, 0, 0},
                                              xfdtd::Vector{0, 16e-3, 40e-3}))};
  constexpr auto max_freq{3e8 / (20 * dz)};
  auto tau{1 / (2 * max_freq)};
  auto t_0{4.5 * tau};
  auto v_s{std::make_shared<xfdtd::VoltageSource>(
      "v_s",
      std::make_unique<xfdtd::Cube>(xfdtd::Vector{-0.787e-3, 0, 24e-3},
                                    xfdtd::Vector{0.787e-3, 0, 2.4e-3}),
      xfdtd::Axis::Direction::XP, 50, xfdtd::Waveform::gaussian(tau, t_0))};

  auto v1{std::make_shared<xfdtd::VoltageMonitor>(
      "v1",
      std::make_unique<xfdtd::Cube>(xfdtd::Vector{-0.787e-3, 0, 24e-3},
                                    xfdtd::Vector{0.787e-3, 0, 2.4e-3}),
      xfdtd::Axis::Direction::XP, "./data/inverted_f_antenna")};

  auto c1{std::make_shared<xfdtd::CurrentMonitor>(
      "c1",
      std::make_unique<xfdtd::Cube>(xfdtd::Vector{-0.39e-3, 0, 24e-3},
                                    xfdtd::Vector{0, 0, 2.4e-3}),
      xfdtd::Axis::Direction::XP, "./data/inverted_f_antenna")};

  auto port_1{std::make_shared<xfdtd::Port>(1, true, 50, c1, v1)};
  auto network{std::make_shared<xfdtd::Network>(
      std::vector<std::shared_ptr<xfdtd::Port>>{port_1},
      xt::arange<double>(2e7, 10e9, 2e7), "./data/inverted_f_antenna")};

  auto nffft{std::make_shared<xfdtd::NFFFT>(13, 13, 13,
                                            xt::xarray<double>{2.4e9, 5.8e9},
                                            "./data/inverted_f_antenna")};

  auto s{xfdtd::Simulation{dx, dy, dz, 0.9}};
  s.addObject(domain);
  s.addObject(substrate);
  s.addObject(antenna_0);
  s.addObject(antenna_1);
  s.addObject(antenna_2);
  s.addObject(antenna_3);
  s.addObject(antenna_4);
  s.addObject(ground);
  s.addObject(v_s);
  s.addMonitor(v1);
  s.addMonitor(c1);
  s.addBoundary(std::make_shared<xfdtd::PML>(8, xfdtd::Axis::Direction::XN));
  s.addBoundary(std::make_shared<xfdtd::PML>(8, xfdtd::Axis::Direction::XP));
  s.addBoundary(std::make_shared<xfdtd::PML>(8, xfdtd::Axis::Direction::YN));
  s.addBoundary(std::make_shared<xfdtd::PML>(8, xfdtd::Axis::Direction::YP));
  s.addBoundary(std::make_shared<xfdtd::PML>(8, xfdtd::Axis::Direction::ZN));
  s.addBoundary(std::make_shared<xfdtd::PML>(8, xfdtd::Axis::Direction::ZP));
  s.addNetwork(network);
  s.addNF2FF(nffft);
  s.run(7000);

  network->output();
  v1->output();
  c1->output();
  nffft->outputRadiationPower();
  nffft->processFarField(
      xfdtd::constant::PI * 0.5,
      xt::linspace<double>(-xfdtd::constant::PI, xfdtd::constant::PI, 360),
      "xy");
  nffft->processFarField(
      xt::linspace<double>(-xfdtd::constant::PI, xfdtd::constant::PI, 360), 0,
      "xz");
  nffft->processFarField(
      xt::linspace<double>(-xfdtd::constant::PI, xfdtd::constant::PI, 360),
      xfdtd::constant::PI * 0.5, "yz");
}

int main() {
  invertedFAntenna();
  return 0;
}
