#include <cstddef>
#include <memory>
#include <xtensor/xnpy.hpp>

#include "xfdtd/boundary/pml.h"
#include "xfdtd/coordinate_system/coordinate_system.h"
#include "xfdtd/material/material.h"
#include "xfdtd/monitor/field_monitor.h"
#include "xfdtd/monitor/movie_monitor.h"
#include "xfdtd/nffft/nffft.h"
#include "xfdtd/object/object.h"
#include "xfdtd/shape/cube.h"
#include "xfdtd/shape/sphere.h"
#include "xfdtd/simulation/simulation.h"
#include "xfdtd/waveform/waveform.h"
#include "xfdtd/waveform_source/tfsf_3d.h"

void dielectricSphereScatter(int num_thread, xfdtd::Divider::Type type) {
  constexpr double dl{5e-3};

  auto domain{std::make_shared<xfdtd::Object>(
      "domain",
      std::make_unique<xfdtd::Cube>(xfdtd::Vector{-0.3, -0.3, -0.3},
                                    xfdtd::Vector{0.6, 0.6, 0.6}),
      xfdtd::Material::createAir())};

  auto dielectric_sphere{std::make_shared<xfdtd::Object>(
      "dielectric_sphere_",
      std::make_unique<xfdtd::Sphere>(xfdtd::Vector{0, 0, 0}, 0.1),
      std::make_unique<xfdtd::Material>(
          "a", xfdtd::ElectroMagneticProperty{3, 2, 0, 0}))};

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
      "./data/dielectric_sphere_scatter")};

  auto movie_xz{std::make_shared<xfdtd::MovieMonitor>(
      std::make_unique<xfdtd::FieldMonitor>(
          std::make_unique<xfdtd::Cube>(xfdtd::Vector{-0.3, -0.3, -0.3},
                                        xfdtd::Vector{0.6, 0.6, 0.6}),
          xfdtd::Axis::XYZ::Y, xfdtd::EMF::Field::EZ, "", ""),
      20, "movie_xz", "./data/dielectric_sphere_scatter/movie_xz")};

  auto movie_yz{std::make_shared<xfdtd::MovieMonitor>(
      std::make_unique<xfdtd::FieldMonitor>(
          std::make_unique<xfdtd::Cube>(xfdtd::Vector{-0.3, -0.3, -0.3},
                                        xfdtd::Vector{0.6, 0.6, 0.6}),
          xfdtd::Axis::XYZ::X, xfdtd::EMF::Field::EZ, "", ""),
      20, "movie_yz", "./data/dielectric_sphere_scatter/movie_yz")};

  auto s{xfdtd::Simulation{dl, dl, dl, 0.9, num_thread, type}};
  s.addObject(domain);
  s.addObject(dielectric_sphere);
  s.addWaveformSource(tfsf);
  s.addNF2FF(nffft_fd);
  s.addBoundary(std::make_shared<xfdtd::PML>(10, xfdtd::Axis::Direction::XN));
  s.addBoundary(std::make_shared<xfdtd::PML>(10, xfdtd::Axis::Direction::XP));
  s.addBoundary(std::make_shared<xfdtd::PML>(10, xfdtd::Axis::Direction::YN));
  s.addBoundary(std::make_shared<xfdtd::PML>(10, xfdtd::Axis::Direction::YP));
  s.addBoundary(std::make_shared<xfdtd::PML>(10, xfdtd::Axis::Direction::ZN));
  s.addBoundary(std::make_shared<xfdtd::PML>(10, xfdtd::Axis::Direction::ZP));
  //   s.addMonitor(movie_xz);
  //   s.addMonitor(movie_yz);
  s.run(2300);

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

  auto time{tfsf->waveform()->time()};
  auto incident_wave_data{tfsf->waveform()->value()};
  xt::dump_npy("./data/dielectric_sphere_scatter/time.npy", time);
  xt::dump_npy("./data/dielectric_sphere_scatter/incident_wave.npy",
               incident_wave_data);
}

int main(int argc, char* argv[]) {
  int num_thread = 1;
  xfdtd::Divider::Type type = xfdtd::Divider::Type::X;
  if (argc > 1) {
    num_thread = std::stoi(argv[1]);
    if (argc > 2) {
      if (std::string(argv[2]) == "Y") {
        type = xfdtd::Divider::Type::Y;
      } else if (std::string(argv[2]) == "Z") {
        type = xfdtd::Divider::Type::Z;
      } else {
        type = xfdtd::Divider::Type::X;
      }
    }
  }
  dielectricSphereScatter(num_thread, type);
  return 0;
}
