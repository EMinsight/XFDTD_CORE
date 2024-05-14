#include <memory>
#include <xtensor/xnpy.hpp>

#include "xfdtd/boundary/pml.h"
#include "xfdtd/common/constant.h"
#include "xfdtd/coordinate_system/coordinate_system.h"
#include "xfdtd/material/material.h"
#include "xfdtd/nffft/nffft_frequency_domain.h"
#include "xfdtd/object/object.h"
#include "xfdtd/shape/cube.h"
#include "xfdtd/simulation/simulation.h"
#include "xfdtd/waveform/waveform.h"
#include "xfdtd/waveform_source/tfsf_3d.h"

void dielectricCubeScatter() {
  constexpr double dl{5e-3};

  auto domain{std::make_shared<xfdtd::Object>(
      "domain",
      std::make_unique<xfdtd::Cube>(
          xfdtd::Vector{-80e-3 - 10 * dl, -80e-3 - 10 * dl, -80e-3 - 10 * dl},
          xfdtd::Vector{160e-3 + 20 * dl, 160e-3 + 20 * dl, 160e-3 + 20 * dl}),
      xfdtd::Material::createAir())};

  auto dielectric_cube{std::make_shared<xfdtd::Object>(
      "dielectric_cube",
      std::make_unique<xfdtd::Cube>(xfdtd::Vector{-80e-3, -80e-3, -80e-3},
                                    xfdtd::Vector{160e-3, 160e-3, 160e-3}),
      std::make_unique<xfdtd::Material>(
          "a", xfdtd::ElectroMagneticProperty{5, 1, 0, 0}))};

  auto l_min{dl * 20};
  auto tau{l_min / 6e8};
  auto t_0{4.5 * tau};
  auto tfsf_3d{std::make_shared<xfdtd::TFSF3D>(
      15, 15, 15, xfdtd::constant::PI * 0.25, xfdtd::constant::PI / 6.0, 0,
      xfdtd::Waveform::gaussian(tau, t_0))};

  auto nffft_fd{std::make_shared<xfdtd::NFFFTFrequencyDomain>(
      13, 13, 13, xt::xarray<double>{1e9},
      "./data/dielectric_cube_scatter/fd")};

  auto s{xfdtd::Simulation{dl, dl, dl, 0.9}};
  s.addObject(domain);
  s.addObject(dielectric_cube);
  s.addWaveformSource(tfsf_3d);
  s.addNF2FF(nffft_fd);
  s.addBoundary(std::make_shared<xfdtd::PML>(10, xfdtd::Axis::Direction::XN));
  s.addBoundary(std::make_shared<xfdtd::PML>(10, xfdtd::Axis::Direction::XP));
  s.addBoundary(std::make_shared<xfdtd::PML>(10, xfdtd::Axis::Direction::YN));
  s.addBoundary(std::make_shared<xfdtd::PML>(10, xfdtd::Axis::Direction::YP));
  s.addBoundary(std::make_shared<xfdtd::PML>(10, xfdtd::Axis::Direction::ZN));
  s.addBoundary(std::make_shared<xfdtd::PML>(10, xfdtd::Axis::Direction::ZP));
  s.run(2000);

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

  auto time{tfsf_3d->waveform()->time()};
  auto incident_wave_data{tfsf_3d->waveform()->value()};
  xt::dump_npy("./data/dielectric_cube_scatter/time.npy", time);
  xt::dump_npy("./data/dielectric_cube_scatter/incident_wave.npy",
               incident_wave_data);
}

int main() {
  dielectricCubeScatter();
  return 0;
}
