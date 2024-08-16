#include <filesystem>

#include "argparse.hpp"
#include "flow_field.h"
#include "xfdtd/boundary/pml.h"
#include "xfdtd/common/type_define.h"
#include "xfdtd/monitor/field_monitor.h"
#include "xfdtd/monitor/movie_monitor.h"
#include "xfdtd/object/object.h"
#include "xfdtd/parallel/parallelized_config.h"
#include "xfdtd/simulation/simulation.h"
#include "xfdtd/waveform_source/tfsf_3d.h"

int main(int argc, char** argv) {
  constexpr auto data_path_str = "./tmp/data/flow_field";
  const auto data_path = std::filesystem::path{data_path_str};

  auto program = argparse::ArgumentParser("flow_field");
  program.add_argument("-f_p", "--flow_field_path")
      .help("flow field data path")
      .required();
  try {
    program.parse_args(argc, argv);
  } catch (const std::runtime_error& err) {
    std::cout << err.what() << std::endl;
    std::cout << program;
    exit(0);
  }

  auto flow_field = std::make_shared<xfdtd::FlowField>(
      "flow_field", program.get<std::string>("--flow_field_path"));

  auto&& shape = flow_field->flowFieldShape();
  auto&& cube = shape.wrappedCube();

  constexpr xfdtd::Real dl{50e-3};

  auto domain_shape = xfdtd::Cube{
      xfdtd::Vector{-10 * dl + cube->originX(), -10 * dl + cube->originY(),
                    -10 * dl + cube->originZ()},
      xfdtd::Vector{20 * dl + cube->sizeX(), 20 * dl + cube->sizeY(),
                    20 * dl + cube->sizeZ()}};

  auto domain = std::make_shared<xfdtd::Object>(
      "domain", std::make_unique<xfdtd::Cube>(domain_shape),
      xfdtd::Material::createAir());

  auto s = xfdtd::Simulation{dl, dl, dl, 0.9, xfdtd::ThreadConfig{4, 1, 1}};
  s.addObject(domain);
  s.addObject(flow_field);

  constexpr auto l_min = dl * 20;
  constexpr auto f_max = 3e8 / l_min;
  constexpr auto tau = l_min / 6e8;
  constexpr auto t_0 = 4.5 * tau;
  constexpr xfdtd::Index tfsf_start = 13;
  auto tfsf{std::make_shared<xfdtd::TFSF3D>(
      tfsf_start, tfsf_start, tfsf_start, xfdtd::constant::PI / 2, 0, 1,
      xfdtd::Waveform::gaussian(tau, t_0))};

  s.addWaveformSource(tfsf);

  auto movie_ez_xy{std::make_shared<xfdtd::MovieMonitor>(
      std::make_unique<xfdtd::FieldMonitor>(
          std::make_unique<xfdtd::Cube>(
              xfdtd::Vector{domain_shape.originX(), domain_shape.originY(), 0},
              xfdtd::Vector{domain_shape.sizeX(), domain_shape.sizeY(), dl}),
          xfdtd::EMF::Field::EZ, "", ""),
      15, "movie_ez_xy", (data_path / "movie_ez_xy").string())};

  s.addMonitor(movie_ez_xy);

  s.addBoundary(std::make_shared<xfdtd::PML>(8, xfdtd::Axis::Direction::XN));
  s.addBoundary(std::make_shared<xfdtd::PML>(8, xfdtd::Axis::Direction::XP));
  s.addBoundary(std::make_shared<xfdtd::PML>(8, xfdtd::Axis::Direction::YN));
  s.addBoundary(std::make_shared<xfdtd::PML>(8, xfdtd::Axis::Direction::YP));
  s.addBoundary(std::make_shared<xfdtd::PML>(8, xfdtd::Axis::Direction::ZN));
  s.addBoundary(std::make_shared<xfdtd::PML>(8, xfdtd::Axis::Direction::ZP));

  s.run(1000);

  return 0;
}