#include <cstddef>
#include <memory>
#include <xtensor/xnpy.hpp>

#include "xfdtd/boundary/pml.h"
#include "xfdtd/coordinate_system/coordinate_system.h"
#include "xfdtd/material/material.h"
#include "xfdtd/monitor/field_monitor.h"
#include "xfdtd/monitor/field_time_monitor.h"
#include "xfdtd/monitor/movie_monitor.h"
#include "xfdtd/nffft/nffft.h"
#include "xfdtd/object/object.h"
#include "xfdtd/parallel/mpi_support.h"
#include "xfdtd/parallel/parallelized_config.h"
#include "xfdtd/shape/cube.h"
#include "xfdtd/shape/sphere.h"
#include "xfdtd/simulation/simulation.h"
#include "xfdtd/waveform/waveform.h"
#include "xfdtd/waveform_source/tfsf_3d.h"

void dielectricSphereScatter() {
  xfdtd::MpiSupport::setMpiParallelDim(2, 2, 1);
  constexpr double dl{7.5e-3};

  auto domain{std::make_shared<xfdtd::Object>(
      "domain",
      std::make_unique<xfdtd::Cube>(xfdtd::Vector{-0.175, -0.175, -0.175},
                                    xfdtd::Vector{0.35, 0.35, 0.35}),
      xfdtd::Material::createAir())};

  auto dielectric_sphere{std::make_shared<xfdtd::Object>(
      "dielectric_sphere",
      std::make_unique<xfdtd::Sphere>(xfdtd::Vector{0, 0, 0}, 0.1),
      std::make_unique<xfdtd::Material>(
          "a", xfdtd::ElectroMagneticProperty{3, 2, 0, 0}))};

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
      nffft_start, nffft_start, nffft_start, xt::xarray<double>{1e9},
      "./data/dielectric_sphere_scatter")};

  auto movie_ex_xz{std::make_shared<xfdtd::MovieMonitor>(
      std::make_unique<xfdtd::FieldMonitor>(
          std::make_unique<xfdtd::Cube>(xfdtd::Vector{-0.175, 0, -0.175},
                                        xfdtd::Vector{0.35, dl, 0.35}),
          xfdtd::Axis::XYZ::Y, xfdtd::EMF::Field::EX, "", ""),
      10, "movie_ex_xz", "./data/dielectric_sphere_scatter/movie_ex_xz")};

  auto movie_ex_yz{std::make_shared<xfdtd::MovieMonitor>(
      std::make_unique<xfdtd::FieldMonitor>(
          std::make_unique<xfdtd::Cube>(xfdtd::Vector{0, -0.175, -0.175},
                                        xfdtd::Vector{dl, 0.35, 0.35}),
          xfdtd::Axis::XYZ::X, xfdtd::EMF::Field::EX),
      10, "movie_ex_yz", "./data/dielectric_sphere_scatter/movie_ex_yz")};

  auto movie_ex_xy = std::make_shared<xfdtd::MovieMonitor>(
      std::make_unique<xfdtd::FieldMonitor>(
          std::make_unique<xfdtd::Cube>(xfdtd::Vector{-0.175, -0.175, 0},
                                        xfdtd::Vector{0.35, 0.35, dl}),
          xfdtd::Axis::XYZ::Z, xfdtd::EMF::Field::EX),
      10, "movie_ex_xy", "./data/dielectric_sphere_scatter/movie_ex_xy");

  auto point_ex = std::make_shared<xfdtd::FieldTimeMonitor>(
      "ex_point",
      std::make_unique<xfdtd::Cube>(xfdtd::Vector{0, -0.175 + 3 * dl, 0},
                                    xfdtd::Vector{dl, dl, dl}),
      "./data/dielectric_sphere_scatter",
      xfdtd::EMF::Field::EX);  // record scatter field

  auto s{xfdtd::Simulation{dl, dl, dl, 0.9, xfdtd::ThreadConfig{2, 2, 1}}};
  s.addObject(domain);
  s.addObject(dielectric_sphere);
  s.addWaveformSource(tfsf);
  s.addNF2FF(nffft_fd);
  s.addBoundary(std::make_shared<xfdtd::PML>(8, xfdtd::Axis::Direction::XN));
  s.addBoundary(std::make_shared<xfdtd::PML>(8, xfdtd::Axis::Direction::XP));
  s.addBoundary(std::make_shared<xfdtd::PML>(8, xfdtd::Axis::Direction::YN));
  s.addBoundary(std::make_shared<xfdtd::PML>(8, xfdtd::Axis::Direction::YP));
  s.addBoundary(std::make_shared<xfdtd::PML>(8, xfdtd::Axis::Direction::ZN));
  s.addBoundary(std::make_shared<xfdtd::PML>(8, xfdtd::Axis::Direction::ZP));
  //   s.addMonitor(movie_ex_xz);
  //   s.addMonitor(movie_ex_yz);
  //   s.addMonitor(movie_ex_xy);
  s.addMonitor(point_ex);
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

  auto time{tfsf->waveform()->time()};
  auto incident_wave_data{tfsf->waveform()->value()};
  xt::dump_npy("./data/dielectric_sphere_scatter/time.npy", time);
  xt::dump_npy("./data/dielectric_sphere_scatter/incident_wave.npy",
               incident_wave_data);
  point_ex->output();
}

int main(int argc, char* argv[]) {
  dielectricSphereScatter();
  return 0;
}
