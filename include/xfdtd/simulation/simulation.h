#ifndef _XFDTD_CORE_SIMULATION_H_
#define _XFDTD_CORE_SIMULATION_H_

#include <xfdtd/boundary/boundary.h>
#include <xfdtd/calculation_param/calculation_param.h>
#include <xfdtd/electromagnetic_field/electromagnetic_field.h>
#include <xfdtd/exception/exception.h>
#include <xfdtd/grid_space/grid_space.h>
#include <xfdtd/monitor/monitor.h>
#include <xfdtd/network/network.h>
#include <xfdtd/nffft/nffft.h>
#include <xfdtd/object/object.h>
#include <xfdtd/parallel/parallelized_config.h>
#include <xfdtd/waveform_source/waveform_source.h>

#include <barrier>
#include <chrono>
#include <cstddef>
#include <memory>
#include <utility>
#include <vector>

namespace xfdtd {

// Forward declaration
class Updator;
class Domain;

class XFDTDSimulationException : public XFDTDException {
 public:
  explicit XFDTDSimulationException(
      std::string message = "XFDTD Simulation Exception")
      : XFDTDException(std::move(message)) {}
};

class Simulation {
 public:
  Simulation(Real dx, Real dy, Real dz, Real cfl,
             ThreadConfig thread_config = ThreadConfig{1, 1, 1});

  ~Simulation();

  bool isRoot() const;

  int myRank() const;

  int numNode() const;

  int numThread() const;

  void addObject(std::shared_ptr<xfdtd::Object> object);

  void addWaveformSource(std::shared_ptr<WaveformSource> waveform_source);

  void addBoundary(std::shared_ptr<Boundary> boundary);

  void addMonitor(std::shared_ptr<Monitor> monitor);

  void addNetwork(std::shared_ptr<Network> network);

  void addNF2FF(std::shared_ptr<NFFFT> nffft);

  void run(std::size_t time_step);

  const std::shared_ptr<CalculationParam>& calculationParam() const;

  const std::shared_ptr<GridSpace>& gridSpace() const;

  const std::shared_ptr<EMF>& emf() const;

  void init();

  auto init(Index time_step) -> void;

  auto waveformSources() -> std::vector<std::shared_ptr<WaveformSource>> {
    return _waveform_sources;
  }

  auto boundaries() -> std::vector<std::shared_ptr<Boundary>> {
    return _boundaries;
  }

  auto monitors() -> std::vector<std::shared_ptr<Monitor>> { return _monitors; }

 private:
  Real _dx, _dy, _dz;
  Real _cfl;
  ThreadConfig _thread_config;
  std::barrier<> _barrier;  // move to thread config

  std::vector<std::shared_ptr<xfdtd::Object>> _objects;
  std::vector<std::shared_ptr<WaveformSource>> _waveform_sources;
  std::vector<std::shared_ptr<Boundary>> _boundaries;
  std::vector<std::shared_ptr<Monitor>> _monitors;
  std::vector<std::shared_ptr<Network>> _networks;
  std::vector<std::shared_ptr<NFFFT>> _nfffts;

  std::chrono::high_resolution_clock::time_point _start_time;
  std::chrono::high_resolution_clock::time_point _end_time;

  std::shared_ptr<GridSpace> _global_grid_space;
  std::shared_ptr<GridSpace> _grid_space;
  std::shared_ptr<CalculationParam> _calculation_param;
  std::shared_ptr<EMF> _emf;

  std::vector<std::unique_ptr<Domain>> _domains;

  std::unique_ptr<CalculationParam> makeCalculationParam();

  std::unique_ptr<TimeParam> makeTimeParam();

  std::unique_ptr<MaterialParam> makeMaterialParam();

  void generateDomain();

  void generateGridSpace();

  void globalGridSpaceDecomposition();

  void generateEMF();

  void generateMaterialSpace();

  void generateFDTDUpdateCoefficient();

  void correctMaterialSpace();

  void correctUpdateCoefficient();

  std::unique_ptr<Updator> makeUpdator(const IndexTask& task);
};

}  // namespace xfdtd

#endif  // _XFDTD_CORE_SIMULATION_H_
