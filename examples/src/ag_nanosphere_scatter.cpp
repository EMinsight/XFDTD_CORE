#include <filesystem>
#include <memory>
#include <xtensor.hpp>
#include <xtensor/xnpy.hpp>

#include "xfdtd/boundary/pml.h"
#include "xfdtd/calculation_param/calculation_param.h"
#include "xfdtd/common/constant.h"
#include "xfdtd/common/type_define.h"
#include "xfdtd/coordinate_system/coordinate_system.h"
#include "xfdtd/material/dispersive_material.h"
#include "xfdtd/monitor/field_monitor.h"
#include "xfdtd/monitor/field_time_monitor.h"
#include "xfdtd/monitor/movie_monitor.h"
#include "xfdtd/nffft/nffft_frequency_domain.h"
#include "xfdtd/nffft/nffft_time_domain.h"
#include "xfdtd/object/object.h"
#include "xfdtd/parallel/mpi_support.h"
#include "xfdtd/shape/cube.h"
#include "xfdtd/shape/sphere.h"
#include "xfdtd/simulation/simulation.h"
#include "xfdtd/waveform_source/tfsf_3d.h"

auto agNanoSphereScatter() {
  xfdtd::MpiSupport::setMpiParallelDim(1, 2, 2);
  constexpr double radius = 40e-9;

  constexpr double dl = 1e-9;

  constexpr auto size = 120 * dl;

  auto&& ag_drude = xfdtd::DrudeMedium::makeDrudeMedium(
      "ag_drude", 3.7, {1.3825e+16}, {2.7347e+13});

  auto&& ag_m_lor = xfdtd::MLorentzMaterial::makeMLorentz(
      "ag_m_lor", 3.7, {1.9114e32}, {0}, {0}, {2.76362e13}, {1});

  const auto domain_area =
      xfdtd::Cube{xfdtd::Vector{-size / 2, -size / 2, -size / 2},
                  xfdtd::Vector{size, size, size}};

  auto domain = std::make_shared<xfdtd::Object>(
      "domain", std::make_unique<xfdtd::Cube>(domain_area),
      xfdtd::Material::createAir());

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

  auto run_simulation = [&](auto&& material) {
    const auto path = std::filesystem::path("./data/ag_nano_sphere_scatter_" +
                                            material->name());

    auto ag_sphere = std::make_shared<xfdtd::Object>(
        "ag_sphere",
        std::make_unique<xfdtd::Sphere>(xfdtd::Vector{0, 0, 0}, radius),
        std::forward<decltype(material)>(material));

    auto ex_point_monitor = std::make_shared<xfdtd::FieldTimeMonitor>(
        "ex_point",
        std::make_unique<xfdtd::Cube>(xfdtd::Vector{0, 0, 0},
                                      xfdtd::Vector{dl, dl, dl}),
        path.string(), xfdtd::EMF::Field::EX);

    auto ex_movie_xz = std::make_shared<xfdtd::MovieMonitor>(
        std::make_unique<xfdtd::FieldMonitor>(
            std::make_unique<xfdtd::Cube>(
                xfdtd::Vector{-size / 2, 0, -size / 2},
                xfdtd::Vector{size, dl, size}),
            xfdtd::EMF::Field::EX, "", ""),
        50, "ex_movie_xz", (path / "movie_xz").string());

    auto nffft_fd = std::make_shared<xfdtd::NFFFTFrequencyDomain>(
        12, 12, 12, xfdtd::Array1D<xfdtd::Real>{800 * 1e12});

    auto nffft_td = std::make_shared<xfdtd::NFFFTTimeDomain>(
        12, 12, 12, xfdtd::constant::PI, 0);

    auto s = xfdtd::Simulation{dl, dl, dl, 0.9, xfdtd::ThreadConfig{1, 1, 1}};
    s.addObject(domain);
    s.addObject(ag_sphere);
    s.addWaveformSource(tfsf_3d);
    s.addMonitor(ex_point_monitor);
    // s.addMonitor(ex_movie_xz);
    s.addNF2FF(nffft_fd);
    s.addNF2FF(nffft_td);
    s.addBoundary(std::make_shared<xfdtd::PML>(8, xfdtd::Axis::Direction::XN));
    s.addBoundary(std::make_shared<xfdtd::PML>(8, xfdtd::Axis::Direction::XP));
    s.addBoundary(std::make_shared<xfdtd::PML>(8, xfdtd::Axis::Direction::YN));
    s.addBoundary(std::make_shared<xfdtd::PML>(8, xfdtd::Axis::Direction::YP));
    s.addBoundary(std::make_shared<xfdtd::PML>(8, xfdtd::Axis::Direction::ZN));
    s.addBoundary(std::make_shared<xfdtd::PML>(8, xfdtd::Axis::Direction::ZP));
    s.run(15000);

    nffft_td->setOutputDir((path / "td").string());
    nffft_td->processFarField();

    nffft_fd->setOutputDir((path / "fd").string());
    nffft_fd->processFarField(
        xt::linspace<xfdtd::Real>(-xfdtd::constant::PI, xfdtd::constant::PI,
                                  360),
        0, "xz");

    if (!s.isRoot()) {
      return;
    }

    auto ag = std::dynamic_pointer_cast<xfdtd::LinearDispersiveMaterial>(
        ag_sphere->material());
    auto epsilon = ag->relativePermittivity(
        xt::linspace<xfdtd::Real>(10 * 1e12, 1200 * 1e12, 100));
    xt::dump_npy(
        (path / "relative_permittivity.npy").string(),
        xt::stack(xt::xtuple(
            xt::linspace<xfdtd::Real>(10 * 1e12, 1200 * 1e12, 100), epsilon)));

    auto time = tfsf_3d->waveform()->time();
    xt::dump_npy((path / "time.npy").string(), time);
    auto incident_wave_data = tfsf_3d->waveform()->value();
    xt::dump_npy((path / "incident_wave.npy").string(),
                 xt::stack(xt::xtuple(time, incident_wave_data)));
    ex_point_monitor->output();
  };

  run_simulation(std::move(ag_drude));

  xfdtd::MpiSupport::instance().barrier();
  if (xfdtd::MpiSupport::instance().isRoot()) {
    std::cout << "Next: ";
  }

  run_simulation(std::move(ag_m_lor));
}

int main() {
  agNanoSphereScatter();
  return 0;
}
