/**
 * @file mini_test.cpp
 * @author your name (you@domain.com)
 * @brief This file is used to test bug in xt::xtuple
 * @version 0.1
 * @date 2024-04-12
 *
 * @copyright Copyright (c) 2024
 *
 */
#include <complex>
#include <memory>
#include <xtensor/xnpy.hpp>

#include "xfdtd/coordinate_system/coordinate_system.h"
#include "xfdtd/material/dispersive_material.h"
#include "xfdtd/material/material.h"
#include "xfdtd/monitor/field_monitor.h"
#include "xfdtd/monitor/movie_monitor.h"
#include "xfdtd/object/object.h"
#include "xfdtd/parallel/parallelized_config.h"
#include "xfdtd/shape/cube.h"
#include "xfdtd/simulation/simulation.h"
#include "xfdtd/waveform/waveform.h"
#include "xfdtd/waveform_source/tfsf_3d.h"

void test() {
  auto dispersive_material =
      std::make_shared<xfdtd::DebyeMedium>(xfdtd::DebyeMedium{
          "debye_medium", 2, {7}, {2e-9 / (2 * xfdtd::constant::PI)}});
  xfdtd::Array1D<double> freq = xt::linspace<double>(0.2e9, 5e9, 100);
  xfdtd::Array1D<std::complex<double>> relative_permittivity =
      dispersive_material->relativePermittivity(freq);
  // Reason: use xt::xtuple
  // Do test in 20 times
  // If add flag -fsanitize, the error will be disappear
  //   xt::dump_npy(
  //       "./data/dispersive_material_scatter/debye_medium/"
  //       "relative_permittivity.npy",
  //       xt::stack(xt::xtuple(freq, relative_permittivity)));

  constexpr double dl{7.5e-3};

  auto domain{std::make_shared<xfdtd::Object>(
      "domain",
      std::make_unique<xfdtd::Cube>(xfdtd::Vector{-0.175, -0.175, -0.175},
                                    xfdtd::Vector{0.35, 0.35, 0.35}),
      xfdtd::Material::createAir())};

  constexpr auto l_min{dl * 20};
  constexpr auto tau{l_min / 6e8};
  constexpr auto t_0{4.5 * tau};
  constexpr std::size_t tfsf_start{static_cast<size_t>(3)};
  auto tfsf{std::make_shared<xfdtd::TFSF3D>(
      tfsf_start, tfsf_start, tfsf_start, 0, 0, 0,
      std::make_unique<xfdtd::Waveform>(xfdtd::Waveform::gaussian(tau, t_0)))};

  auto movie_ex_xz{std::make_shared<xfdtd::MovieMonitor>(
      std::make_unique<xfdtd::FieldMonitor>(
          std::make_unique<xfdtd::Cube>(xfdtd::Vector{-0.175, 0, -0.175},
                                        xfdtd::Vector{0.35, dl, 0.35}),
          xfdtd::Axis::XYZ::Y, xfdtd::EMF::Field::EX, "", ""),
      10, "movie_ex_xz",
      "./data/dispersive_material_scatter/debye_medium/movie_ex_xz")};

  auto s{xfdtd::Simulation{dl, dl, dl, 0.9, xfdtd::ThreadConfig{1, 1, 1}}};
  s.addObject(domain);
  s.addWaveformSource(tfsf);
  s.addMonitor(movie_ex_xz);
  s.run(120);
}

int main(int argc, char* argv[]) {
  test();
  return 0;
}
