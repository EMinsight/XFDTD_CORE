#include <filesystem>
#include <memory>
#include <xtensor/xnpy.hpp>

#include "xfdtd/coordinate_system/coordinate_system.h"
#include "xfdtd/material/material.h"
#include "xfdtd/monitor/voltage_monitor.h"
#include "xfdtd/object/lumped_element/capacitor.h"
#include "xfdtd/object/lumped_element/pec_plane.h"
#include "xfdtd/object/lumped_element/voltage_source.h"
#include "xfdtd/object/object.h"
#include "xfdtd/parallel/parallelized_config.h"
#include "xfdtd/shape/cube.h"
#include "xfdtd/simulation/simulation.h"
#include "xfdtd/waveform/waveform.h"

constexpr static xfdtd::Real SIZE{1e-3};

constexpr static xfdtd::Real ORIGIN_A{-1 * SIZE};
constexpr static xfdtd::Real ORIGIN_B{-1 * SIZE};
constexpr static xfdtd::Real ORIGIN_C{-1 * SIZE};
constexpr static xfdtd::Real LENGTH_A{6 * SIZE};
constexpr static xfdtd::Real LENGTH_B{4 * SIZE};
constexpr static xfdtd::Real LENGTH_C{4 * SIZE};

constexpr static xfdtd::Real PLANE_O_A{0 * SIZE};
constexpr static xfdtd::Real PLANE_O_B{0 * SIZE};
constexpr static xfdtd::Real PLANE_O_C{0 * SIZE};
constexpr static xfdtd::Real PLANE_L_A{1 * SIZE};
constexpr static xfdtd::Real PLANE_L_B{1 * SIZE};
constexpr static xfdtd::Real PLANE_L_C{0 * SIZE};

constexpr static xfdtd::Real PLANE1_O_A{0 * SIZE};
constexpr static xfdtd::Real PLANE1_O_B{0 * SIZE};
constexpr static xfdtd::Real PLANE1_O_C{1 * SIZE};
constexpr static xfdtd::Real PLANE1_L_A{1 * SIZE};
constexpr static xfdtd::Real PLANE1_L_B{1 * SIZE};
constexpr static xfdtd::Real PLANE1_L_C{0 * SIZE};

constexpr static xfdtd::Real V_SOURCE_O_A{0 * SIZE};
constexpr static xfdtd::Real V_SOURCE_O_B{0 * SIZE};
constexpr static xfdtd::Real V_SOURCE_O_C{0 * SIZE};
constexpr static xfdtd::Real V_SOURCE_L_A{0 * SIZE};
constexpr static xfdtd::Real V_SOURCE_L_B{1 * SIZE};
constexpr static xfdtd::Real V_SOURCE_L_C{1 * SIZE};

constexpr static xfdtd::Real CAPACITOR_O_A{1 * SIZE};
constexpr static xfdtd::Real CAPACITOR_O_B{0 * SIZE};
constexpr static xfdtd::Real CAPACITOR_O_C{0 * SIZE};
constexpr static xfdtd::Real CAPACITOR_L_A{0 * SIZE};
constexpr static xfdtd::Real CAPACITOR_L_B{1 * SIZE};
constexpr static xfdtd::Real CAPACITOR_L_C{1 * SIZE};

constexpr static xfdtd::Real V_MONITOR_O_A{1 * SIZE};
constexpr static xfdtd::Real V_MONITOR_O_B{0 * SIZE};
constexpr static xfdtd::Real V_MONITOR_O_C{0 * SIZE};
constexpr static xfdtd::Real V_MONITOR_L_A{0 * SIZE};
constexpr static xfdtd::Real V_MONITOR_L_B{1 * SIZE};
constexpr static xfdtd::Real V_MONITOR_L_C{1 * SIZE};

const static std::string OUT_DIR{"./tmp/data/simple_capacitor/"};

xfdtd::Real exactSimpleCapacitorSolver(xfdtd::Real c, xfdtd::Real r,
                                       xfdtd::Real t, xfdtd::Real t_0) {
  t = t - t_0;
  if (t < 0) {
    return 0;
  }
  return 1 - std::exp(-t / (r * c));
}

void capacitorX() {
  // X: C A B
  auto domain{std::make_shared<xfdtd::Object>(
      "domain",
      std::make_unique<xfdtd::Cube>(
          xfdtd::Vector{ORIGIN_C, ORIGIN_A, ORIGIN_B},
          xfdtd::Vector{LENGTH_C, LENGTH_A, LENGTH_B}),
      xfdtd::Material::createAir())};

  auto plane{std::make_shared<xfdtd::PecPlane>(
      "xn_plane", std::make_unique<xfdtd::Cube>(
                      xfdtd::Vector{PLANE_O_C, PLANE_O_A, PLANE_O_B},
                      xfdtd::Vector{PLANE_L_C, PLANE_L_A, PLANE_L_B}))};

  auto plane2{std::make_shared<xfdtd::PecPlane>(
      "xp_plane", std::make_unique<xfdtd::Cube>(
                      xfdtd::Vector{PLANE1_O_C, PLANE1_O_A, PLANE1_O_B},
                      xfdtd::Vector{PLANE1_L_C, PLANE1_L_A, PLANE1_L_B}))};

  auto dt{xfdtd::TimeParam::calculateDt(0.98, SIZE, SIZE, SIZE)};
  auto v_source{std::make_shared<xfdtd::VoltageSource>(
      "v",
      std::make_unique<xfdtd::Cube>(
          xfdtd::Vector{V_SOURCE_O_C, V_SOURCE_O_A, V_SOURCE_O_B},
          xfdtd::Vector{V_SOURCE_L_C, V_SOURCE_L_A, V_SOURCE_L_B}),
      xfdtd::Axis::Direction::XP, 50, xfdtd::Waveform::step(100 * dt))};

  auto capacitor{std::make_shared<xfdtd::Capacitor>(
      "capacitor",
      std::make_unique<xfdtd::Cube>(
          xfdtd::Vector{CAPACITOR_O_C, CAPACITOR_O_A, CAPACITOR_O_B},
          xfdtd::Vector{CAPACITOR_L_C, CAPACITOR_L_A, CAPACITOR_L_B}),
      xfdtd::Axis::XYZ::X, 1e-11)};

  auto v_monitor{std::make_shared<xfdtd::VoltageMonitor>(
      "v_monitor",
      std::make_unique<xfdtd::Cube>(
          xfdtd::Vector{V_MONITOR_O_C, V_MONITOR_O_A, V_MONITOR_O_B},
          xfdtd::Vector{V_MONITOR_L_C, V_MONITOR_L_A, V_MONITOR_L_B}),
      xfdtd::Axis::Direction::XP, OUT_DIR)};

  auto simulation{xfdtd::Simulation{SIZE, SIZE, SIZE, 0.98}};
  simulation.addObject(domain);
  simulation.addObject(plane);
  simulation.addObject(plane2);
  simulation.addObject(v_source);
  simulation.addObject(capacitor);
  simulation.addMonitor(v_monitor);
  simulation.run(2000);

  v_monitor->output();

  auto h_time{simulation.calculationParam()->timeParam()->hTime()};
  auto exact{xt::zeros_like(h_time)};
  for (std::size_t i{0}; i < h_time.size(); ++i) {
    exact(i) = exactSimpleCapacitorSolver(1e-11, 50, h_time(i), 100 * dt);
  }
  xt::dump_npy((std::filesystem::path{OUT_DIR} / "exact.npy").string(), exact);
}

void capacitorZ() {
  auto domain{std::make_shared<xfdtd::Object>(
      "domain",
      std::make_unique<xfdtd::Cube>(
          xfdtd::Vector{ORIGIN_A, ORIGIN_B, ORIGIN_C},
          xfdtd::Vector{LENGTH_A, LENGTH_B, LENGTH_C}),
      xfdtd::Material::createAir())};

  auto plane{std::make_shared<xfdtd::PecPlane>(
      "zn_plane", std::make_unique<xfdtd::Cube>(
                      xfdtd::Vector{PLANE_O_A, PLANE_O_B, PLANE_O_C},
                      xfdtd::Vector{PLANE_L_A, PLANE_L_B, PLANE_L_C}))};

  auto plane2{std::make_shared<xfdtd::PecPlane>(
      "zp_plane", std::make_unique<xfdtd::Cube>(
                      xfdtd::Vector{PLANE1_O_A, PLANE1_O_B, PLANE1_O_C},
                      xfdtd::Vector{PLANE1_L_A, PLANE1_L_B, PLANE1_L_C}))};
  auto dt{xfdtd::TimeParam::calculateDt(0.98, SIZE, SIZE, SIZE)};
  auto v_source{std::make_shared<xfdtd::VoltageSource>(
      "v",
      std::make_unique<xfdtd::Cube>(
          xfdtd::Vector{V_SOURCE_O_A, V_SOURCE_O_B, V_SOURCE_O_C},
          xfdtd::Vector{V_SOURCE_L_A, V_SOURCE_L_B, V_SOURCE_L_C}),
      xfdtd::Axis::Direction::ZP, 50, xfdtd::Waveform::step(100 * dt))};

  auto capacitor{std::make_shared<xfdtd::Capacitor>(
      "capacitor",
      std::make_unique<xfdtd::Cube>(
          xfdtd::Vector{CAPACITOR_O_A, CAPACITOR_O_B, CAPACITOR_O_C},
          xfdtd::Vector{CAPACITOR_L_A, CAPACITOR_L_B, CAPACITOR_L_C}),
      xfdtd::Axis::XYZ::Z, 1e-11)};

  auto v_monitor{std::make_shared<xfdtd::VoltageMonitor>(
      "v_monitor",
      std::make_unique<xfdtd::Cube>(
          xfdtd::Vector{V_MONITOR_O_A, V_MONITOR_O_B, V_MONITOR_O_C},
          xfdtd::Vector{V_MONITOR_L_A, V_MONITOR_L_B, V_MONITOR_L_C}),
      xfdtd::Axis::Direction::ZP, OUT_DIR)};

  auto simulation{
      xfdtd::Simulation{SIZE, SIZE, SIZE, 0.98, xfdtd::ThreadConfig{1, 1, 1}}};
  simulation.addObject(domain);
  simulation.addObject(plane);
  simulation.addObject(plane2);
  simulation.addObject(v_source);
  simulation.addObject(capacitor);
  simulation.addMonitor(v_monitor);
  simulation.run(2000);

  v_monitor->output();

  auto h_time{simulation.calculationParam()->timeParam()->hTime()};
  auto exact{xt::zeros_like(h_time)};
  for (std::size_t i{0}; i < h_time.size(); ++i) {
    exact(i) = exactSimpleCapacitorSolver(1e-11, 50, h_time(i), 100 * dt);
  }

  xt::dump_npy((std::filesystem::path{OUT_DIR} / "exact.npy").string(), exact);
}

int main(int argc, char* argv[]) {
  capacitorX();
  return 0;
}
