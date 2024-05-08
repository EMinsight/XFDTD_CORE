#include <complex>
#include <filesystem>
#include <memory>
#include <xtensor.hpp>
#include <xtensor/xnpy.hpp>

#include "xfdtd/boundary/pml.h"
#include "xfdtd/calculation_param/calculation_param.h"
#include "xfdtd/common/type_define.h"
#include "xfdtd/coordinate_system/coordinate_system.h"
#include "xfdtd/material/dispersive_material.h"
#include "xfdtd/monitor/field_monitor.h"
#include "xfdtd/monitor/field_time_monitor.h"
#include "xfdtd/monitor/movie_monitor.h"
#include "xfdtd/object/object.h"
#include "xfdtd/parallel/mpi_support.h"
#include "xfdtd/shape/cube.h"
#include "xfdtd/shape/sphere.h"
#include "xfdtd/simulation/simulation.h"
#include "xfdtd/waveform/waveform.h"
#include "xfdtd/waveform_source/tfsf_3d.h"

auto nanoSphereScatter() -> void {
  xfdtd::MpiSupport::setMpiParallelDim(1, 2, 2);
  const auto path = std::filesystem::path("./data/au_nano_sphere_scatter");

  constexpr double radius = 40e-9;

  constexpr double dl = 1e-9;

  constexpr auto size = 150 * dl;

  const auto domain_area =
      xfdtd::Cube{xfdtd::Vector{-size / 2, -size / 2, -size / 2},
                  xfdtd::Vector{size, size, size}};

  auto domain = std::make_shared<xfdtd::Object>(
      "domain", std::make_unique<xfdtd::Cube>(domain_area),
      xfdtd::Material::createAir());

  using Complex = std::complex<xfdtd::Real>;
  auto p = xfdtd::Array1D<Complex>{
      Complex{-5.6502e+05, -1.4392e+14}, Complex{-2.6071e+05, -1.5755e+15},
      Complex{-3.2315e+14, -3.8103e+15}, Complex{-1.6773e+15, -3.5824e+15}};
  auto r = xfdtd::Array1D<Complex>{
      Complex{3.5065e+15, 5.6365e+17}, Complex{2.0571e+14, 1.0699e+14},
      Complex{4.7948e+14, -3.2664e+13}, Complex{9.3298e+15, -4.4758e+14}};

  auto&& a_0 = -2 * xt::real(r * xt::conj(p));
  auto&& a_1 = 2 * xt::real(r);
  auto&& b_0 = xt::pow(xt::abs(p), 2);
  auto&& b_1 = -2 * xt::real(p);
  auto&& b_2 = xt::ones_like(b_1);

  auto au_sphere = std::make_shared<xfdtd::Object>(
      "au_sphere",
      std::make_unique<xfdtd::Sphere>(xfdtd::Vector{0, 0, 0}, radius),
      xfdtd::MLorentzMaterial::makeMLorentz("Au_CCPR", 1.366, a_0, a_1, b_0,
                                            b_1, b_2));

  constexpr auto l_min = dl / 20;
  constexpr auto f = 3e8 / l_min;
  constexpr auto f_max = 800e12 * 1.2;
  constexpr auto f_min = 200e12 * 1.2;
  constexpr auto bandwidth = (f_max - f_min);
  constexpr auto tau = 0.996 / bandwidth;
  constexpr auto t_0 = 4.5 * tau;
  auto&& waveform =
      xfdtd::Waveform::cosineModulatedGaussian(tau, t_0, (f_max + f_min) / 2);
  auto tfsf_3d =
      std::make_shared<xfdtd::TFSF3D>(15, 15, 15, 0, 0, 0, std::move(waveform));

  auto ex_point_monitor = std::make_shared<xfdtd::FieldTimeMonitor>(
      "ex_point",
      std::make_unique<xfdtd::Cube>(xfdtd::Vector{0, 0, 0},
                                    xfdtd::Vector{dl, dl, dl}),
      path.string(), xfdtd::EMF::Field::EX);

  auto ex_movie_yz = std::make_shared<xfdtd::MovieMonitor>(
      std::make_unique<xfdtd::FieldMonitor>(
          std::make_unique<xfdtd::Cube>(xfdtd::Vector{0, -size / 2, -size / 2},
                                        xfdtd::Vector{dl, size, size}),
          xfdtd::EMF::Field::EX, "", ""),
      10, "ex_movie_yz", (path / "movie_yz").string());

  auto s = xfdtd::Simulation{dl, dl, dl, 0.9, xfdtd::ThreadConfig{2, 1, 1}};
  s.addObject(domain);
  s.addObject(au_sphere);
  s.addWaveformSource(tfsf_3d);
  s.addMonitor(ex_point_monitor);
  s.addMonitor(ex_movie_yz);
  s.addBoundary(std::make_shared<xfdtd::PML>(8, xfdtd::Axis::Direction::XN));
  s.addBoundary(std::make_shared<xfdtd::PML>(8, xfdtd::Axis::Direction::XP));
  s.addBoundary(std::make_shared<xfdtd::PML>(8, xfdtd::Axis::Direction::YN));
  s.addBoundary(std::make_shared<xfdtd::PML>(8, xfdtd::Axis::Direction::YP));
  s.addBoundary(std::make_shared<xfdtd::PML>(8, xfdtd::Axis::Direction::ZN));
  s.addBoundary(std::make_shared<xfdtd::PML>(8, xfdtd::Axis::Direction::ZP));
  s.run(15000);

  if (!s.isRoot()) {
    return;
  }

  auto dt = s.calculationParam()->timeParam()->dt();
  auto fft_freq = xt::linspace<double>(0, 1 / (2 * dt), 1000);
  auto au = std::dynamic_pointer_cast<xfdtd::LinearDispersiveMaterial>(
      au_sphere->material());
  auto epsilon = au->relativePermittivity(fft_freq);
  xt::dump_npy((path / "relative_permittivity.npy").string(),
               xt::stack(xt::xtuple(fft_freq, epsilon)));

  auto time = tfsf_3d->waveform()->time();
  auto incident_wave_data = tfsf_3d->waveform()->value();
  xt::dump_npy((path / "incident_wave.npy").string(),
               xt::stack(xt::xtuple(time, incident_wave_data)));
  ex_point_monitor->output();
}

int main(int argc, char** argv) {
  nanoSphereScatter();
  return 0;
}
