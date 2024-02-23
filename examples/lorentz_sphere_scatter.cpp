#include <xfdtd/material/dispersive_material.h>

#include <filesystem>
#include <memory>
#include <xtensor/xnpy.hpp>

#include "xfdtd/boundary/pml.h"
#include "xfdtd/coordinate_system/coordinate_system.h"
#include "xfdtd/object/object.h"
#include "xfdtd/shape/cube.h"
#include "xfdtd/shape/sphere.h"
#include "xfdtd/simulation/simulation.h"
#include "xfdtd/util/constant.h"
#include "xfdtd/waveform_source/tfsf_3d.h"

static const std::filesystem::path DATA_DIR = "./data";

void lorentzSphereScatter() {
  const auto lorentz_sphere_scatter_dir = DATA_DIR / "lorentz_sphere_scatter";
  if (!std::filesystem::exists(lorentz_sphere_scatter_dir)) {
    std::filesystem::create_directories(lorentz_sphere_scatter_dir);
  }

  constexpr double dl{5e-3};

  auto domain = std::make_shared<xfdtd::Object>(
      "domain",
      std::make_unique<xfdtd::Cube>(xfdtd::Vector{-0.3, -0.3, -0.3},
                                    xfdtd::Vector{0.6, 0.6, 0.6}),
      xfdtd::Material::createAir());

  auto sphere_material = std::make_shared<xfdtd::LorentzMedium>(
      xfdtd::LorentzMedium{"scatter_m",
                           2,
                           {5},
                           {2 * xfdtd::constant::PI * 2e9},
                           {xfdtd::constant::PI * 2e9}});

  constexpr double radius = 1e-1;
  auto sphere = std::make_shared<xfdtd::Object>(
      "scatter",
      std::make_unique<xfdtd::Sphere>(xfdtd::Vector{0, 0, 0}, radius),
      sphere_material);

  constexpr auto l_min{dl * 20};
  constexpr auto tau{l_min / 6e8};
  constexpr auto t_0{4.5 * tau};
  constexpr std::size_t left_len{static_cast<size_t>((0.3 / dl))};
  constexpr std::size_t remain_len{static_cast<size_t>(left_len - 0.1 / dl)};
  constexpr std::size_t tfsf_start{static_cast<size_t>((remain_len / 2) + 10)};
  auto tfsf{std::make_shared<xfdtd::TFSF3D>(
      tfsf_start, tfsf_start, tfsf_start, 0, 0, 0,
      std::make_unique<xfdtd::Waveform>(xfdtd::Waveform::gaussian(tau, t_0)))};

  constexpr std::size_t nffft_start{static_cast<size_t>((remain_len / 3) + 10)};
  auto nffft_fd{std::make_shared<xfdtd::NFFFT>(
      nffft_start, nffft_start, nffft_start, xt::xarray<double>{1e9},
      (lorentz_sphere_scatter_dir).string())};

  auto simulation = xfdtd::Simulation{dl, dl, dl, 0.99};
  simulation.addObject(domain);
  simulation.addObject(sphere);
  simulation.addWaveformSource(tfsf);
  simulation.addNF2FF(nffft_fd);
  simulation.addBoundary(
      std::make_shared<xfdtd::PML>(10, xfdtd::Axis::Direction::XN));
  simulation.addBoundary(
      std::make_shared<xfdtd::PML>(10, xfdtd::Axis::Direction::XP));
  simulation.addBoundary(
      std::make_shared<xfdtd::PML>(10, xfdtd::Axis::Direction::YN));
  simulation.addBoundary(
      std::make_shared<xfdtd::PML>(10, xfdtd::Axis::Direction::YP));
  simulation.addBoundary(
      std::make_shared<xfdtd::PML>(10, xfdtd::Axis::Direction::ZN));
  simulation.addBoundary(
      std::make_shared<xfdtd::PML>(10, xfdtd::Axis::Direction::ZP));
  simulation.run(3000);

  auto relative_permittivity =
      sphere_material->relativePermittivity(xt::linspace<double>(0, 5e9, 1000));

  nffft_fd->processFarField(
      xfdtd::constant::PI * 0.5,
      xt::linspace<double>(-xfdtd::constant::PI, xfdtd::constant::PI, 360),
      "xy");
  xt::dump_npy(
      (lorentz_sphere_scatter_dir / "relative_permittivity.npy").string(),
      relative_permittivity);
  nffft_fd->processFarField(
      xt::linspace<double>(-xfdtd::constant::PI, xfdtd::constant::PI, 360), 0,
      "xz");
  nffft_fd->processFarField(
      xt::linspace<double>(-xfdtd::constant::PI, xfdtd::constant::PI, 360),
      xfdtd::constant::PI * 0.5, "yz");

  auto time = tfsf->waveform()->time();
  auto incident_wave_data = tfsf->waveform()->value();
  xt::dump_npy((lorentz_sphere_scatter_dir / "time.npy").string(), time);
  xt::dump_npy((lorentz_sphere_scatter_dir / "incident_wave.npy").string(),
               incident_wave_data);
}

int main() {
  lorentzSphereScatter();
  return 0;
}
