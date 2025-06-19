#include <memory>
#include <xtensor/xnpy.hpp>

#include "xfdtd/coordinate_system/coordinate_system.h"
#include "xfdtd/material/material.h"
#include "xfdtd/monitor/voltage_monitor.h"
#include "xfdtd/object/lumped_element/capacitor.h"
#include "xfdtd/object/lumped_element/inductor.h"
#include "xfdtd/object/lumped_element/pec_plane.h"
#include "xfdtd/object/lumped_element/voltage_source.h"
#include "xfdtd/object/object.h"
#include "xfdtd/parallel/mpi_support.h"
#include "xfdtd/parallel/parallelized_config.h"
#include "xfdtd/shape/cube.h"
#include "xfdtd/simulation/simulation.h"
#include "xfdtd/waveform/waveform.h"

constexpr xfdtd::Real SIZE{1e-3};

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
constexpr static xfdtd::Real PLANE1_O_C{2 * SIZE};

constexpr static xfdtd::Real PLANE1_L_A{1 * SIZE};
constexpr static xfdtd::Real PLANE1_L_B{1 * SIZE};
constexpr static xfdtd::Real PLANE1_L_C{0 * SIZE};

constexpr static xfdtd::Real V_SOURCE_O_A{0 * SIZE};
constexpr static xfdtd::Real V_SOURCE_O_B{0 * SIZE};
constexpr static xfdtd::Real V_SOURCE_O_C{0 * SIZE};
constexpr static xfdtd::Real V_SOURCE_L_A{0 * SIZE};
constexpr static xfdtd::Real V_SOURCE_L_B{1 * SIZE};
constexpr static xfdtd::Real V_SOURCE_L_C{2 * SIZE};

constexpr static xfdtd::Real INDUCTOR_O_A{1 * SIZE};
constexpr static xfdtd::Real INDUCTOR_O_B{0 * SIZE};
constexpr static xfdtd::Real INDUCTOR_O_C{0 * SIZE};
constexpr static xfdtd::Real INDUCTOR_L_A{0 * SIZE};
constexpr static xfdtd::Real INDUCTOR_L_B{1 * SIZE};
constexpr static xfdtd::Real INDUCTOR_L_C{1 * SIZE};

constexpr static xfdtd::Real CAPACITOR_O_A{1 * SIZE};
constexpr static xfdtd::Real CAPACITOR_O_B{0 * SIZE};
constexpr static xfdtd::Real CAPACITOR_O_C{1 * SIZE};
constexpr static xfdtd::Real CAPACITOR_L_A{0 * SIZE};
constexpr static xfdtd::Real CAPACITOR_L_B{1 * SIZE};
constexpr static xfdtd::Real CAPACITOR_L_C{1 * SIZE};

constexpr static xfdtd::Real V_MONITOR_O_A{1 * SIZE};
constexpr static xfdtd::Real V_MONITOR_O_B{0 * SIZE};
constexpr static xfdtd::Real V_MONITOR_O_C{0 * SIZE};
constexpr static xfdtd::Real V_MONITOR_L_A{0 * SIZE};
constexpr static xfdtd::Real V_MONITOR_L_B{1 * SIZE};
constexpr static xfdtd::Real V_MONITOR_L_C{2 * SIZE};

void rlcCircuit() {
  // Y : B C A
  auto domain{std::make_shared<xfdtd::Object>(
      "domain",
      std::make_unique<xfdtd::Cube>(
          xfdtd::Vector{ORIGIN_A, ORIGIN_B, ORIGIN_C},
          xfdtd::Vector{LENGTH_A, LENGTH_B, LENGTH_C}),
      xfdtd::Material::createAir())};

  auto plane{std::make_shared<xfdtd::PecPlane>(
      "plane", std::make_unique<xfdtd::Cube>(
                   xfdtd::Vector{PLANE_O_A, PLANE_O_B, PLANE_O_C},
                   xfdtd::Vector{PLANE_L_A, PLANE_L_B, PLANE_L_C}))};
  auto plane1{std::make_shared<xfdtd::PecPlane>(
      "plane1", std::make_unique<xfdtd::Cube>(
                    xfdtd::Vector{PLANE1_O_A, PLANE1_O_B, PLANE1_O_C},
                    xfdtd::Vector{PLANE1_L_A, PLANE1_L_B, PLANE1_L_C}))};

  constexpr xfdtd::Real bandwidth{4e9};
  constexpr xfdtd::Real tau{0.996 / bandwidth};
  constexpr xfdtd::Real t_0{4.5 * tau};
  auto v_source{std::make_shared<xfdtd::VoltageSource>(
      "v_source",
      std::make_unique<xfdtd::Cube>(
          xfdtd::Vector{V_SOURCE_O_A, V_SOURCE_O_B, V_SOURCE_O_C},
          xfdtd::Vector{V_SOURCE_L_A, V_SOURCE_L_B, V_SOURCE_L_C}),
      xfdtd::Axis::Direction::ZP, 50,
      xfdtd::Waveform::cosineModulatedGaussian(tau, t_0, 2e9))};

  auto inductor{std::make_shared<xfdtd::Inductor>(
      "inductor",
      std::make_unique<xfdtd::Cube>(
          xfdtd::Vector{INDUCTOR_O_A, INDUCTOR_O_B, INDUCTOR_O_C},
          xfdtd::Vector{INDUCTOR_L_A, INDUCTOR_L_B, INDUCTOR_L_C}),
      xfdtd::Axis::XYZ::Z, 1e-8)};

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
      xfdtd::Axis::Direction::ZP, "./tmp/data/rlc_circuit")};

  auto s{
      xfdtd::Simulation{SIZE, SIZE, SIZE, 0.98, xfdtd::ThreadConfig{2, 1, 1}}};
  s.addObject(domain);
  s.addObject(plane);
  s.addObject(plane1);
  s.addObject(v_source);
  s.addObject(inductor);
  s.addObject(capacitor);
  s.addMonitor(v_monitor);
  s.run(2000);

  v_monitor->output();
  xt::dump_npy("./tmp/data/rlc_circuit/source.npy",
               v_source->waveform()->value());
  xt::dump_npy("./tmp/data/rlc_circuit/time.npy",
               s.calculationParam()->timeParam()->hTime());
}

int main(int argc, char *argv[]) {
  xfdtd::MpiSupport::instance(argc, argv);
  rlcCircuit();
}
