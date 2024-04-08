#ifndef __XFDTD_EXAMPLE_DISPERSIVE_SPHERE_SCATTER_H__
#define __XFDTD_EXAMPLE_DISPERSIVE_SPHERE_SCATTER_H__

#include <xfdtd/boundary/pml.h>
#include <xfdtd/common/constant.h>
#include <xfdtd/coordinate_system/coordinate_system.h>
#include <xfdtd/material/dispersive_material.h>
#include <xfdtd/material/material.h>
#include <xfdtd/object/object.h>
#include <xfdtd/shape/cube.h>
#include <xfdtd/shape/sphere.h>
#include <xfdtd/simulation/simulation.h>
#include <xfdtd/waveform_source/tfsf_3d.h>

#include <complex>
#include <filesystem>
#include <memory>
#include <string>
#include <xtensor/xnpy.hpp>

inline constexpr double FREQ = 1e9;

inline auto relativePermittivityToNonDispersive(
    const std::complex<double>& relative_permittivity, const double& freq) {
  const auto epsilon_r = std::real(relative_permittivity);
  const auto sigma = std::abs(std::imag(relative_permittivity)) *
                     xfdtd::constant::EPSILON_0 * 2 * xfdtd::constant::PI *
                     freq;
  return xfdtd::ElectroMagneticProperty{epsilon_r, 1, sigma, 0};
}

inline auto getNonDispersiveMaterial(
    double freq, std::complex<double> relative_permittivity) {
  auto em = relativePermittivityToNonDispersive(relative_permittivity, freq);
  return std::make_unique<xfdtd::Material>("non_dispersive", em);
}

inline void outputRelativePermittivity(
    const xt::xarray<double>& freq,
    const std::shared_ptr<xfdtd::LinearDispersiveMaterial>& material,
    const std::string& file_name) {
  xt::dump_npy(file_name, material->relativePermittivity(freq));
}

inline void runSimulation(std::shared_ptr<xfdtd::Material> sphere_material,
                          std::string_view dir) {
  const std::filesystem::path sphere_scatter_dir{dir};
  if (!std::filesystem::exists(sphere_scatter_dir) ||
      !std::filesystem::is_directory(sphere_scatter_dir)) {
    std::filesystem::create_directories(sphere_scatter_dir);
  }

  std::cout << "Save dir: " << std::filesystem::absolute(sphere_scatter_dir)
            << "\n";

  if (sphere_material->dispersion()) {
    outputRelativePermittivity(
        xt::linspace<double>(0, 5 * FREQ, 1000),
        std::dynamic_pointer_cast<xfdtd::LinearDispersiveMaterial>(
            sphere_material),
        (sphere_scatter_dir / "relative_permittivity.npy").string());
  }

  constexpr double dl{7.5e-3};

  auto domain{std::make_shared<xfdtd::Object>(
      "domain",
      std::make_unique<xfdtd::Cube>(xfdtd::Vector{-0.175, -0.175, -0.175},
                                    xfdtd::Vector{0.35, 0.35, 0.35}),
      xfdtd::Material::createAir())};

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
  constexpr std::size_t tfsf_start{static_cast<size_t>(15)};
  auto tfsf{std::make_shared<xfdtd::TFSF3D>(
      tfsf_start, tfsf_start, tfsf_start, 0, 0, 0,
      std::make_unique<xfdtd::Waveform>(xfdtd::Waveform::gaussian(tau, t_0)))};

  constexpr std::size_t nffft_start{static_cast<size_t>(11)};
  auto nffft_fd{std::make_shared<xfdtd::NFFFT>(
      nffft_start, nffft_start, nffft_start, xt::xarray<double>{FREQ},
      (sphere_scatter_dir).string())};

  auto simulation = xfdtd::Simulation{dl, dl, dl, 0.9};
  simulation.addObject(domain);
  simulation.addObject(sphere);
  simulation.addWaveformSource(tfsf);
  simulation.addNF2FF(nffft_fd);
  simulation.addBoundary(
      std::make_shared<xfdtd::PML>(8, xfdtd::Axis::Direction::XN));
  simulation.addBoundary(
      std::make_shared<xfdtd::PML>(8, xfdtd::Axis::Direction::XP));
  simulation.addBoundary(
      std::make_shared<xfdtd::PML>(8, xfdtd::Axis::Direction::YN));
  simulation.addBoundary(
      std::make_shared<xfdtd::PML>(8, xfdtd::Axis::Direction::YP));
  simulation.addBoundary(
      std::make_shared<xfdtd::PML>(8, xfdtd::Axis::Direction::ZN));
  simulation.addBoundary(
      std::make_shared<xfdtd::PML>(8, xfdtd::Axis::Direction::ZP));

  simulation.run(2200);

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

  auto time = tfsf->waveform()->time();
  auto incident_wave_data = tfsf->waveform()->value();
  xt::dump_npy((sphere_scatter_dir / "time.npy").string(), time);
  xt::dump_npy((sphere_scatter_dir / "incident_wave.npy").string(),
               incident_wave_data);
}

inline void testCase(
    const std::shared_ptr<xfdtd::LinearDispersiveMaterial>& material, int id) {
  std::filesystem::path data_dir =
      std::filesystem::path("./data/dispersive_material_scatter");

  if (id != 0) {
    std::cout << "Dispersive material: " + material->name() + "\n";
    data_dir /= material->name();

    runSimulation(material, data_dir.string());
  } else {
    auto non_dispersive_material = getNonDispersiveMaterial(
        FREQ, material->relativePermittivity({FREQ}).front());
    std::cout << "Non-dispersive material: " + material->name() + "\n";

    data_dir /= "non_dispersive_" + material->name();
    std::cout << "ElectroMagneticProperty: "
              << non_dispersive_material->emProperty().toString() << "\n";
    runSimulation(std::move(non_dispersive_material), data_dir.string());
  }
}

#endif  // __XFDTD_EXAMPLE_DISPERSIVE_SPHERE_SCATTER_H__
