#include <iostream>
#include <memory>
#include <sstream>
#include <vector>

#include "argparse.hpp"
#include "xfdtd/boundary/pml.h"
#include "xfdtd/common/type_define.h"
#include "xfdtd/coordinate_system/coordinate_system.h"
#include "xfdtd/material/material.h"
#include "xfdtd/object/object.h"
#include "xfdtd/parallel/mpi_support.h"
#include "xfdtd/parallel/parallelized_config.h"
#include "xfdtd/shape/cube.h"
#include "xfdtd/simulation/simulation.h"

void dielectricSphereScatter(int argc, char* argv[]) {
  auto program = argparse::ArgumentParser("benchmark");
  program.add_argument("-t", "--time_steps")
      .help("Number of time steps")
      .default_value(4800)
      .scan<'d', int>();

  program.add_argument("-d", "--dl")
      .help("Grid resolution")
      .default_value(xfdtd::Real{5e-3})
      .scan<'g', xfdtd::Real>();

  program.add_argument("-m_c", "--mpi_config")
      .help("MPI configuration")
      .default_value(std::vector<int>{2, 2, 1})
      .nargs(3)
      .scan<'d', int>();

  program.add_argument("-t_c", "--thread_config")
      .help("Thread configuration")
      .default_value(std::vector<int>{1, 1, 1})
      .nargs(3)
      .scan<'d', int>();

  program.add_argument("--with_pml")
      .help("With PML")
      .default_value(false)
      .implicit_value(true);

  try {
    program.parse_args(argc, argv);
  } catch (const std::runtime_error& err) {
    std::stringstream ss;
    ss << err.what() << '\n';
    ss << program << '\n';
    std::cerr << ss.str();
    xfdtd::MpiSupport::instance().abort(1);
  }

  auto time_steps = program.get<int>("-t");
  auto dl = program.get<xfdtd::Real>("-d");
  auto mpi_config = program.get<std::vector<int>>("-m_c");
  auto thread_config = program.get<std::vector<int>>("-t_c");
  auto with_pml = program.get<bool>("--with_pml");

  xfdtd::MpiSupport::setMpiParallelDim(mpi_config[0], mpi_config[1],
                                       mpi_config[2]);
  if (xfdtd::MpiSupport::instance().isRoot()) {
    std::cout << "time_steps: " << time_steps << '\n';
    std::cout << "dl: " << dl << '\n';
    std::cout << "mpi_config: " << mpi_config[0] << " " << mpi_config[1] << " "
              << mpi_config[2] << '\n';
    std::cout << "thread_config: " << thread_config[0] << " "
              << thread_config[1] << " " << thread_config[2] << '\n';
    std::cout << "with_pml: " << std::boolalpha << with_pml << '\n';
  }

  auto domain{std::make_shared<xfdtd::Object>(
      "domain",
      std::make_unique<xfdtd::Cube>(xfdtd::Vector{-0.175, -0.175, -0.175},
                                    xfdtd::Vector{0.35, 0.35, 0.35}),
      xfdtd::Material::createAir())};

  auto s{
      xfdtd::Simulation{dl, dl, dl, 0.9,
                        xfdtd::ThreadConfig{thread_config[0], thread_config[1],
                                            thread_config[2]}}};
  s.addObject(domain);
  if (with_pml) {
    s.addBoundary(std::make_shared<xfdtd::PML>(8, xfdtd::Axis::Direction::XN));
    s.addBoundary(std::make_shared<xfdtd::PML>(8, xfdtd::Axis::Direction::XP));
    s.addBoundary(std::make_shared<xfdtd::PML>(8, xfdtd::Axis::Direction::YN));
    s.addBoundary(std::make_shared<xfdtd::PML>(8, xfdtd::Axis::Direction::YP));
    s.addBoundary(std::make_shared<xfdtd::PML>(8, xfdtd::Axis::Direction::ZN));
    s.addBoundary(std::make_shared<xfdtd::PML>(8, xfdtd::Axis::Direction::ZP));
  }
  s.addDefaultVisitor();
  s.run(time_steps);
}

int main(int argc, char* argv[]) {
  xfdtd::MpiSupport::init(argc, argv);
  dielectricSphereScatter(argc, argv);
  return 0;
}
