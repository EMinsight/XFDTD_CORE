#include <cstddef>
#include <filesystem>
#include <memory>
#include <string>
#include <xtensor/xnpy.hpp>

#include "xfdtd/boundary/pml.h"
#include "xfdtd/coordinate_system/coordinate_system.h"
#include "xfdtd/material/material.h"
#include "xfdtd/object/object.h"
#include "xfdtd/parallel/mpi_support.h"
#include "xfdtd/parallel/parallelized_config.h"
#include "xfdtd/shape/cube.h"
#include "xfdtd/shape/sphere.h"
#include "xfdtd/simulation/simulation.h"
#include "xfdtd/waveform/waveform.h"
#include "xfdtd/waveform_source/tfsf_3d.h"

void dielectricSphereScatter(xfdtd::ThreadConfig thread_config) {
  constexpr double dl{2.5e-3};
  using namespace std::string_view_literals;
  constexpr auto data_path_str = "./data/dielectric_sphere_scatter"sv;
  const auto data_path = std::filesystem::path{data_path_str};

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
  constexpr auto f_max{3e8 / l_min};  // max frequency: 5 GHz in dl = 3e-3
  constexpr auto tau{l_min / 6e8};
  constexpr auto t_0{4.5 * tau};
  constexpr std::size_t tfsf_start{static_cast<size_t>(15)};
  auto tfsf{
      std::make_shared<xfdtd::TFSF3D>(tfsf_start, tfsf_start, tfsf_start, 0, 0,
                                      0, xfdtd::Waveform::gaussian(tau, t_0))};

  auto s{xfdtd::Simulation{dl, dl, dl, 0.9, thread_config}};
  s.addObject(domain);
  s.addObject(dielectric_sphere);
  s.addWaveformSource(tfsf);

  s.addBoundary(std::make_shared<xfdtd::PML>(8, xfdtd::Axis::Direction::XN));
  s.addBoundary(std::make_shared<xfdtd::PML>(8, xfdtd::Axis::Direction::XP));
  s.addBoundary(std::make_shared<xfdtd::PML>(8, xfdtd::Axis::Direction::YN));
  s.addBoundary(std::make_shared<xfdtd::PML>(8, xfdtd::Axis::Direction::YP));
  s.addBoundary(std::make_shared<xfdtd::PML>(8, xfdtd::Axis::Direction::ZN));
  s.addBoundary(std::make_shared<xfdtd::PML>(8, xfdtd::Axis::Direction::ZP));

  s.run(400);
}

int main(int argc, char* argv[]) {
  xfdtd::MpiSupport::setMpiParallelDim(2, 2, 1);
  xfdtd::MpiSupport::instance(argc, argv);

  auto num_x = 1;
  auto num_y = 1;
  auto num_z = 1;
  if (4 <= argc) {
    num_x = std::stoi(argv[1]);
    num_y = std::stoi(argv[2]);
    num_z = std::stoi(argv[3]);
  }
  auto thread_config = xfdtd::ThreadConfig{num_x, num_y, num_z};

  dielectricSphereScatter(thread_config);
  return 0;
}
