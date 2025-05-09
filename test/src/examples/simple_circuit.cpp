#include <memory>

#include "xfdtd/coordinate_system/coordinate_system.h"
#include "xfdtd/material/material.h"
#include "xfdtd/monitor/current_monitor.h"
#include "xfdtd/monitor/voltage_monitor.h"
#include "xfdtd/object/lumped_element/current_source.h"
#include "xfdtd/object/lumped_element/pec_plane.h"
#include "xfdtd/object/lumped_element/resistor.h"
#include "xfdtd/object/lumped_element/voltage_source.h"
#include "xfdtd/object/object.h"
#include "xfdtd/parallel/parallelized_config.h"
#include "xfdtd/shape/cube.h"
#include "xfdtd/simulation/simulation.h"
#include "xfdtd/waveform/waveform.h"

const static std::string OUT_DIR{"./tmp/data/simple_circuit/"};

constexpr static xfdtd::Real SIZE{1e-3};

constexpr static xfdtd::Real ORIGIN_A{-3 * SIZE};
constexpr static xfdtd::Real ORIGIN_B{-3 * SIZE};
constexpr static xfdtd::Real ORIGIN_C{-3 * SIZE};

constexpr static xfdtd::Real LENGTH_A{14 * SIZE};
constexpr static xfdtd::Real LENGTH_B{8 * SIZE};
constexpr static xfdtd::Real LENGTH_C{10 * SIZE};

constexpr static xfdtd::Real PLANE_O_A{0 * SIZE};
constexpr static xfdtd::Real PLANE_O_B{0 * SIZE};
constexpr static xfdtd::Real PLANE_O_C{0 * SIZE};

constexpr static xfdtd::Real PLANE_L_A{8 * SIZE};
constexpr static xfdtd::Real PLANE_L_B{2 * SIZE};
constexpr static xfdtd::Real PLANE_L_C{0 * SIZE};

constexpr static xfdtd::Real PLANE1_O_A{0 * SIZE};
constexpr static xfdtd::Real PLANE1_O_B{0 * SIZE};
constexpr static xfdtd::Real PLANE1_O_C{4 * SIZE};

constexpr static xfdtd::Real PLANE1_L_A{8 * SIZE};
constexpr static xfdtd::Real PLANE1_L_B{2 * SIZE};
constexpr static xfdtd::Real PLANE1_L_C{0 * SIZE};

constexpr static xfdtd::Real V_SOURCE_O_A{0 * SIZE};
constexpr static xfdtd::Real V_SOURCE_O_B{0 * SIZE};
constexpr static xfdtd::Real V_SOURCE_O_C{0 * SIZE};

constexpr static xfdtd::Real V_SOURCE_L_A{0 * SIZE};
constexpr static xfdtd::Real V_SOURCE_L_B{2 * SIZE};
constexpr static xfdtd::Real V_SOURCE_L_C{4 * SIZE};

constexpr static xfdtd::Real RESISTOR_O_A{8 * SIZE};
constexpr static xfdtd::Real RESISTOR_O_B{0 * SIZE};
constexpr static xfdtd::Real RESISTOR_O_C{0 * SIZE};

constexpr static xfdtd::Real RESISTOR_L_A{0 * SIZE};
constexpr static xfdtd::Real RESISTOR_L_B{2 * SIZE};
constexpr static xfdtd::Real RESISTOR_L_C{4 * SIZE};

constexpr static xfdtd::Real V_MONITOR_O_A{5 * SIZE};
constexpr static xfdtd::Real V_MONITOR_O_B{0 * SIZE};
constexpr static xfdtd::Real V_MONITOR_O_C{0 * SIZE};

constexpr static xfdtd::Real V_MONITOR_L_A{0 * SIZE};
constexpr static xfdtd::Real V_MONITOR_L_B{2 * SIZE};
constexpr static xfdtd::Real V_MONITOR_L_C{4 * SIZE};

constexpr static xfdtd::Real I_MONITOR_O_A{5 * SIZE};
constexpr static xfdtd::Real I_MONITOR_O_B{0 * SIZE};
constexpr static xfdtd::Real I_MONITOR_O_C{3 * SIZE};

constexpr static xfdtd::Real I_MONITOR_L_A{0 * SIZE};
constexpr static xfdtd::Real I_MONITOR_L_B{2 * SIZE};
constexpr static xfdtd::Real I_MONITOR_L_C{1 * SIZE};

void simpleCircuitX() {
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

  auto v_source{std::make_shared<xfdtd::CurrentSource>(
      "v",
      std::make_unique<xfdtd::Cube>(
          xfdtd::Vector{V_SOURCE_O_C, V_SOURCE_O_A, V_SOURCE_O_B},
          xfdtd::Vector{V_SOURCE_L_C, V_SOURCE_L_A, V_SOURCE_L_B}),
      xfdtd::Axis::Direction::XP, 0, xfdtd::Waveform::sine(5e8))};

  auto resistor{std::make_shared<xfdtd::Resistor>(
      "r",
      std::make_unique<xfdtd::Cube>(
          xfdtd::Vector{RESISTOR_O_C, RESISTOR_O_A, RESISTOR_O_B},
          xfdtd::Vector{RESISTOR_L_C, RESISTOR_L_A, RESISTOR_L_B}),
      xfdtd::Axis::XYZ::X, 50)};

  auto circuit_v_monitor{std::make_shared<xfdtd::VoltageMonitor>(
      "v_monitor",
      std::make_unique<xfdtd::Cube>(
          xfdtd::Vector{V_MONITOR_O_C, V_MONITOR_O_A, V_MONITOR_O_B},
          xfdtd::Vector{V_MONITOR_L_C, V_MONITOR_L_A, V_MONITOR_L_B}),
      xfdtd::Axis::Direction::XP, OUT_DIR)};

  auto circuit_i_monitor{std::make_shared<xfdtd::CurrentMonitor>(
      "i_monitor",
      std::make_unique<xfdtd::Cube>(
          xfdtd::Vector{I_MONITOR_O_C, I_MONITOR_O_A, I_MONITOR_O_B},
          xfdtd::Vector{I_MONITOR_L_C, I_MONITOR_L_A, I_MONITOR_L_B}),
      xfdtd::Axis::Direction::YP, OUT_DIR)};

  auto s{xfdtd::Simulation{SIZE, SIZE, SIZE, 0.98}};
  s.addObject(domain);
  s.addObject(plane);
  s.addObject(plane2);
  s.addObject(v_source);
  s.addObject(resistor);
  s.addMonitor(circuit_v_monitor);
  s.addMonitor(circuit_i_monitor);
  s.run(3000);

  circuit_i_monitor->output();
  circuit_v_monitor->output();
}

void simpleCircuitY() {
  auto domain{std::make_shared<xfdtd::Object>(
      "domain",
      std::make_unique<xfdtd::Cube>(
          xfdtd::Vector{ORIGIN_B, ORIGIN_C, ORIGIN_A},
          xfdtd::Vector{LENGTH_B, LENGTH_C, LENGTH_A}),
      xfdtd::Material::createAir())};

  auto plane{std::make_shared<xfdtd::PecPlane>(
      "yn_plane", std::make_unique<xfdtd::Cube>(
                      xfdtd::Vector{PLANE_O_B, PLANE_O_C, PLANE_O_A},
                      xfdtd::Vector{PLANE_L_B, PLANE_L_C, PLANE_L_A}))};

  auto plane2{std::make_shared<xfdtd::PecPlane>(
      "yp_plane", std::make_unique<xfdtd::Cube>(
                      xfdtd::Vector{PLANE1_O_B, PLANE1_O_C, PLANE1_O_A},
                      xfdtd::Vector{PLANE1_L_B, PLANE1_L_C, PLANE1_L_A}))};

  auto v_source{std::make_shared<xfdtd::VoltageSource>(
      "v",
      std::make_unique<xfdtd::Cube>(
          xfdtd::Vector{V_SOURCE_O_B, V_SOURCE_O_C, V_SOURCE_O_A},
          xfdtd::Vector{V_SOURCE_L_B, V_SOURCE_L_C, V_SOURCE_L_A}),
      xfdtd::Axis::Direction::YP, 50, xfdtd::Waveform::sine(5e8))};

  auto resistor{std::make_shared<xfdtd::Resistor>(
      "r",
      std::make_unique<xfdtd::Cube>(
          xfdtd::Vector{RESISTOR_O_B, RESISTOR_O_C, RESISTOR_O_A},
          xfdtd::Vector{RESISTOR_L_B, RESISTOR_L_C, RESISTOR_L_A}),
      xfdtd::Axis::XYZ::Y, 150)};

  auto circuit_v_monitor{std::make_shared<xfdtd::VoltageMonitor>(
      "v_monitor",
      std::make_unique<xfdtd::Cube>(
          xfdtd::Vector{V_MONITOR_O_B, V_MONITOR_O_C, V_MONITOR_O_A},
          xfdtd::Vector{V_MONITOR_L_B, V_MONITOR_L_C, V_MONITOR_L_A}),
      xfdtd::Axis::Direction::YP, OUT_DIR)};

  auto circuit_i_monitor{std::make_shared<xfdtd::CurrentMonitor>(
      "i_monitor",
      std::make_unique<xfdtd::Cube>(
          xfdtd::Vector{I_MONITOR_O_B, I_MONITOR_O_C, I_MONITOR_O_A},
          xfdtd::Vector{I_MONITOR_L_B, I_MONITOR_L_C, I_MONITOR_L_A}),
      xfdtd::Axis::Direction::ZP, OUT_DIR)};

  auto s{xfdtd::Simulation{SIZE, SIZE, SIZE, 0.9}};
  s.addObject(domain);
  s.addObject(plane);
  s.addObject(plane2);
  s.addObject(v_source);
  s.addObject(resistor);
  s.addMonitor(circuit_v_monitor);
  s.addMonitor(circuit_i_monitor);
  s.run(3000);

  circuit_i_monitor->output();
  circuit_v_monitor->output();
}

void simpleCircuitZ() {
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

  auto v_source{std::make_shared<xfdtd::VoltageSource>(
      "v",
      std::make_unique<xfdtd::Cube>(
          xfdtd::Vector{V_SOURCE_O_A, V_SOURCE_O_B, V_SOURCE_O_C},
          xfdtd::Vector{V_SOURCE_L_A, V_SOURCE_L_B, V_SOURCE_L_C}),
      xfdtd::Axis::Direction::ZP, 50, xfdtd::Waveform::sine(5e8))};

  auto resistor{std::make_shared<xfdtd::Resistor>(
      "r",
      std::make_unique<xfdtd::Cube>(
          xfdtd::Vector{RESISTOR_O_A, RESISTOR_O_B, RESISTOR_O_C},
          xfdtd::Vector{RESISTOR_L_A, RESISTOR_L_B, RESISTOR_L_C}),
      xfdtd::Axis::XYZ::Z, 50)};

  auto circuit_v_monitor{std::make_shared<xfdtd::VoltageMonitor>(
      "v_monitor",
      std::make_unique<xfdtd::Cube>(
          xfdtd::Vector{V_MONITOR_O_A, V_MONITOR_O_B, V_MONITOR_O_C},
          xfdtd::Vector{V_MONITOR_L_A, V_MONITOR_L_B, V_MONITOR_L_C}),
      xfdtd::Axis::Direction::ZP, OUT_DIR)};

  auto circuit_i_monitor{std::make_shared<xfdtd::CurrentMonitor>(
      "i_monitor",
      std::make_unique<xfdtd::Cube>(
          xfdtd::Vector{I_MONITOR_O_A, I_MONITOR_O_B, I_MONITOR_O_C},
          xfdtd::Vector{I_MONITOR_L_A, I_MONITOR_L_B, I_MONITOR_L_C}),
      xfdtd::Axis::Direction::XP, OUT_DIR)};

  auto s{
      xfdtd::Simulation{SIZE, SIZE, SIZE, 0.90, xfdtd::ThreadConfig{1, 1, 1}}};

  s.addObject(domain);
  s.addObject(plane);
  s.addObject(plane2);
  s.addObject(v_source);
  s.addObject(resistor);
  s.addMonitor(circuit_v_monitor);
  s.addMonitor(circuit_i_monitor);
  s.run(5000);

  circuit_i_monitor->output();
  circuit_v_monitor->output();
}

int main(int argc, char* argv[]) { simpleCircuitZ(); }
