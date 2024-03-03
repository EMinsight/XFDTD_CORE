#include <memory>
#include <xtensor/xnpy.hpp>

#include "divider/divider.h"
#include "xfdtd/coordinate_system/coordinate_system.h"
#include "xfdtd/material/material.h"
#include "xfdtd/monitor/current_monitor.h"
#include "xfdtd/monitor/voltage_monitor.h"
#include "xfdtd/object/lumped_element/current_source.h"
#include "xfdtd/object/lumped_element/pec_plane.h"
#include "xfdtd/object/lumped_element/resistor.h"
#include "xfdtd/object/lumped_element/voltage_source.h"
#include "xfdtd/object/object.h"
#include "xfdtd/shape/cube.h"
#include "xfdtd/simulation/simulation.h"
#include "xfdtd/waveform/waveform.h"

constexpr static double SIZE{1e-3};

constexpr static double ORIGIN_A{-3 * SIZE};
constexpr static double ORIGIN_B{-3 * SIZE};
constexpr static double ORIGIN_C{-3 * SIZE};

constexpr static double LENGTH_A{14 * SIZE};
constexpr static double LENGTH_B{8 * SIZE};
constexpr static double LENGTH_C{10 * SIZE};

constexpr static double PLANE_O_A{0 * SIZE};
constexpr static double PLANE_O_B{0 * SIZE};
constexpr static double PLANE_O_C{0 * SIZE};

constexpr static double PLANE_L_A{8 * SIZE};
constexpr static double PLANE_L_B{2 * SIZE};
constexpr static double PLANE_L_C{0 * SIZE};

constexpr static double PLANE1_O_A{0 * SIZE};
constexpr static double PLANE1_O_B{0 * SIZE};
constexpr static double PLANE1_O_C{4 * SIZE};

constexpr static double PLANE1_L_A{8 * SIZE};
constexpr static double PLANE1_L_B{2 * SIZE};
constexpr static double PLANE1_L_C{0 * SIZE};

constexpr static double V_SOURCE_O_A{0 * SIZE};
constexpr static double V_SOURCE_O_B{0 * SIZE};
constexpr static double V_SOURCE_O_C{0 * SIZE};

constexpr static double V_SOURCE_L_A{0 * SIZE};
constexpr static double V_SOURCE_L_B{2 * SIZE};
constexpr static double V_SOURCE_L_C{4 * SIZE};

constexpr static double RESISTOR_O_A{8 * SIZE};
constexpr static double RESISTOR_O_B{0 * SIZE};
constexpr static double RESISTOR_O_C{0 * SIZE};

constexpr static double RESISTOR_L_A{0 * SIZE};
constexpr static double RESISTOR_L_B{2 * SIZE};
constexpr static double RESISTOR_L_C{4 * SIZE};

constexpr static double V_MONITOR_O_A{5 * SIZE};
constexpr static double V_MONITOR_O_B{0 * SIZE};
constexpr static double V_MONITOR_O_C{0 * SIZE};

constexpr static double V_MONITOR_L_A{0 * SIZE};
constexpr static double V_MONITOR_L_B{2 * SIZE};
constexpr static double V_MONITOR_L_C{4 * SIZE};

constexpr static double I_MONITOR_O_A{5 * SIZE};
constexpr static double I_MONITOR_O_B{0 * SIZE};
constexpr static double I_MONITOR_O_C{3 * SIZE};

constexpr static double I_MONITOR_L_A{0 * SIZE};
constexpr static double I_MONITOR_L_B{2 * SIZE};
constexpr static double I_MONITOR_L_C{1 * SIZE};

void simpleCircuitX(int num_thread) {
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
      xfdtd::Axis::Direction::XP, 0,
      std::make_unique<xfdtd::Waveform>(xfdtd::Waveform::sine(5e8)))};

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
      xfdtd::Axis::Direction::XP, "./data/simple_circuit")};

  auto circuit_i_monitor{std::make_shared<xfdtd::CurrentMonitor>(
      "i_monitor",
      std::make_unique<xfdtd::Cube>(
          xfdtd::Vector{I_MONITOR_O_C, I_MONITOR_O_A, I_MONITOR_O_B},
          xfdtd::Vector{I_MONITOR_L_C, I_MONITOR_L_A, I_MONITOR_L_B}),
      xfdtd::Axis::Direction::YP, "./data/simple_circuit")};

  auto s{xfdtd::Simulation{SIZE, SIZE, SIZE, 0.98, num_thread,
                           xfdtd::Divider::Type::X}};
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

void simpleCircuitY(int num_thread) {
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
      xfdtd::Axis::Direction::YP, 50,
      std::make_unique<xfdtd::Waveform>(xfdtd::Waveform::sine(5e8)))};

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
      xfdtd::Axis::Direction::YP, "./data/simple_circuit")};

  auto circuit_i_monitor{std::make_shared<xfdtd::CurrentMonitor>(
      "i_monitor",
      std::make_unique<xfdtd::Cube>(
          xfdtd::Vector{I_MONITOR_O_B, I_MONITOR_O_C, I_MONITOR_O_A},
          xfdtd::Vector{I_MONITOR_L_B, I_MONITOR_L_C, I_MONITOR_L_A}),
      xfdtd::Axis::Direction::ZP, "./data/simple_circuit")};

  auto s{xfdtd::Simulation{SIZE, SIZE, SIZE, 0.9, num_thread,
                           xfdtd::Divider::Type::Y}};
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

void simpleCircuitZ(int num_thread) {
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
      xfdtd::Axis::Direction::ZP, 50,
      std::make_unique<xfdtd::Waveform>(xfdtd::Waveform::sine(5e8)))};

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
      xfdtd::Axis::Direction::ZP, "./data/simple_circuit")};

  auto circuit_i_monitor{std::make_shared<xfdtd::CurrentMonitor>(
      "i_monitor",
      std::make_unique<xfdtd::Cube>(
          xfdtd::Vector{I_MONITOR_O_A, I_MONITOR_O_B, I_MONITOR_O_C},
          xfdtd::Vector{I_MONITOR_L_A, I_MONITOR_L_B, I_MONITOR_L_C}),
      xfdtd::Axis::Direction::XP, "./data/simple_circuit")};

  auto s{xfdtd::Simulation{SIZE, SIZE, SIZE, 0.90, num_thread,
                           xfdtd::Divider::Type::Z}};

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

int main(int argc, char* argv[]) {
  int num_thread = 1;
  if (argc > 1) {
    if (argc == 3) {
      num_thread = std::stoi(argv[2]);
    }
    if (std::string(argv[1]) == "X") {
      simpleCircuitX(num_thread);
    } else if (std::string(argv[1]) == "Y") {
      simpleCircuitY(num_thread);
    } else if (std::string(argv[1]) == "Z") {
      simpleCircuitZ(num_thread);
    }

  } else {
    simpleCircuitZ(num_thread);
  }
}
