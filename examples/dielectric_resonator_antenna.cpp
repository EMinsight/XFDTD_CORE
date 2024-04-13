#include <memory>
#include <utility>

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
#include "xfdtd/object/lumped_element/voltage_source.h"
#include "xfdtd/object/object.h"
#include "xfdtd/shape/cube.h"
#include "xfdtd/simulation/simulation.h"
#include "xfdtd/waveform/waveform.h"

void dielectricResonatorAntenna() {
  constexpr double dx{0.715e-3};
  constexpr double dy{0.508e-3};
  constexpr double dz{0.5e-3};

  auto domain_cube{std::make_unique<xfdtd::Cube>(
      xfdtd::Vector{-8e-3 - 10 * dx, -6e-3 - 10 * dy, -10 * dz},
      xfdtd::Vector{30e-3 + 20 * dx, 37e-3 + 20 * dy, 26.1e-3 + 20 * dz})};
  auto domain{std::make_shared<xfdtd::Object>(
      "domain",
      std::make_unique<xfdtd::Cube>(
          xfdtd::Vector{-8e-3 - 10 * dx, -6e-3 - 10 * dy, -10 * dz},
          xfdtd::Vector{30e-3 + 20 * dx, 37e-3 + 20 * dy, 26.1e-3 + 20 * dz}),
      xfdtd::Material::createAir())};

  auto box{std::make_shared<xfdtd::Object>(
      "box",
      std::make_unique<xfdtd::Cube>(xfdtd::Vector{0, 0, 0},
                                    xfdtd::Vector{14.3e-3, 25.4e-3, 26.1e-3}),
      std::make_unique<xfdtd::Material>(
          "box", xfdtd::ElectroMagneticProperty{9.8, 1, 0, 0}))};

  auto ground{std::make_shared<xfdtd::PecPlane>(
      "ground", std::make_unique<xfdtd::Cube>(xfdtd::Vector{-8e-3, -6e-3, 0},
                                              xfdtd::Vector{30e-3, 37e-3, 0}))};

  auto strip{std::make_shared<xfdtd::PecPlane>(
      "strip", std::make_unique<xfdtd::Cube>(xfdtd::Vector{0, 12.2e-3, 1e-3},
                                             xfdtd::Vector{0, 1e-3, 9e-3}))};

  auto l_min{dz * 30};
  auto tau{l_min / 6e8};
  auto t_0{4.5 * tau};
  auto v_s{std::make_shared<xfdtd::VoltageSource>(
      "v_s",
      std::make_unique<xfdtd::Cube>(xfdtd::Vector{0, 12.2e-3, 0},
                                    xfdtd::Vector{0, 1e-3, 1e-3}),
      xfdtd::Axis::Direction::ZP, 50,
      std::make_unique<xfdtd::Waveform>(xfdtd::Waveform::gaussian(tau, t_0)))};

  auto v1{std::make_shared<xfdtd::VoltageMonitor>(
      "v1",
      std::make_unique<xfdtd::Cube>(xfdtd::Vector{0, 12.2e-3, 0},
                                    xfdtd::Vector{0, 1e-3, 1e-3}),
      xfdtd::Axis::Direction::ZP, "")};

  auto c1{std::make_shared<xfdtd::CurrentMonitor>(
      "c1",
      std::make_unique<xfdtd::Cube>(xfdtd::Vector{0, 12.2e-3, 0.5e-3},
                                    xfdtd::Vector{0, 1e-3, 0}),
      xfdtd::Axis::Direction::ZP, "")};

  auto e_monitor{std::make_unique<xfdtd::FieldMonitor>(
      std::make_unique<xfdtd::Cube>(
          xfdtd::Vector{0, -6e-3 - 10 * dy, -10 * dz},
          xfdtd::Vector{dx, 37e-3 + 20 * dy, 26.1e-3 + 20 * dz}),
      xfdtd::EMF::Field::EZ, "e", "")};

  auto movie{std::make_shared<xfdtd::MovieMonitor>(
      std::move(e_monitor), 10, "movie",
      "./data/dielectric_resonator_antenna/movie")};

  auto port_1{std::make_shared<xfdtd::Port>(1, true, 50, c1, v1)};

  auto network{std::make_shared<xfdtd::Network>(
      std::vector<std::shared_ptr<xfdtd::Port>>{port_1},
      xt::linspace(2e9, 6e9, 100), "./data/dielectric_resonator_antenna")};

  auto nffft{std::make_shared<xfdtd::NFFFT>(
      13, 13, 13, xt::xarray<double>{3.5e9, 4.3e9},
      "./data/dielectric_resonator_antenna")};

  auto s{xfdtd::Simulation{dx, dy, dz, 0.9}};
  s.addObject(domain);
  s.addObject(box);
  s.addObject(ground);
  s.addObject(strip);
  s.addObject(v_s);
  s.addMonitor(v1);
  s.addMonitor(c1);
  //   s.addMonitor(movie);
  s.addNetwork(network);
  s.addNF2FF(nffft);
  s.addBoundary(std::make_shared<xfdtd::PML>(8, xfdtd::Axis::Direction::XN));
  s.addBoundary(std::make_shared<xfdtd::PML>(8, xfdtd::Axis::Direction::XP));
  s.addBoundary(std::make_shared<xfdtd::PML>(8, xfdtd::Axis::Direction::YN));
  s.addBoundary(std::make_shared<xfdtd::PML>(8, xfdtd::Axis::Direction::YP));
  s.addBoundary(std::make_shared<xfdtd::PML>(8, xfdtd::Axis::Direction::ZN));
  s.addBoundary(std::make_shared<xfdtd::PML>(8, xfdtd::Axis::Direction::ZP));
  s.run(5000);

  network->output();
  v1->setOutputDir("./data/dielectric_resonator_antenna");
  v1->output();
  c1->setOutputDir("./data/dielectric_resonator_antenna");
  c1->output();
  nffft->outputRadiationPower();
  nffft->processFarField(
      xfdtd::constant::PI * 0.5,
      xt::linspace<double>(-xfdtd::constant::PI, xfdtd::constant::PI, 360),
      "xy_plane", domain_cube->center());
  nffft->processFarField(
      xt::linspace<double>(-xfdtd::constant::PI, xfdtd::constant::PI, 360), 0,
      "xz_plane", domain_cube->center());
  nffft->processFarField(
      xt::linspace<double>(-xfdtd::constant::PI, xfdtd::constant::PI, 360),
      xfdtd::constant::PI * 0.5, "yz_plane", domain_cube->center());
}

int main(int argc, char* argv[]) {
  dielectricResonatorAntenna();
  return 0;
}
