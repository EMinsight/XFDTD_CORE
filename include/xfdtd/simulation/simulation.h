#ifndef _XFDTD_LIB_SIMULATION_H_
#define _XFDTD_LIB_SIMULATION_H_

#include <xfdtd/grid_space/grid_space.h>

#include <barrier>
#include <chrono>
#include <cstddef>
#include <memory>
#include <utility>
#include <vector>

#include "divider/divider.h"
#include "domain/domain.h"
#include "xfdtd/boundary/boundary.h"
#include "xfdtd/calculation_param/calculation_param.h"
#include "xfdtd/electromagnetic_field/electromagnetic_field.h"
#include "xfdtd/exception/exception.h"
#include "xfdtd/monitor/monitor.h"
#include "xfdtd/network/network.h"
#include "xfdtd/nffft/nffft.h"
#include "xfdtd/object/object.h"
#include "xfdtd/waveform_source/waveform_source.h"

namespace xfdtd {

// Forward declaration
class Updator;

class XFDTDSimulationException : public XFDTDException {
 public:
  explicit XFDTDSimulationException(
      std::string message = "XFDTD Simulation Exception")
      : XFDTDException(std::move(message)) {}
};

class Simulation {
 public:
  Simulation(double dx, double dy, double dz, double cfl, int num_thread = 1,
             Divider::Type divider_type = Divider::Type::UNDEFINED);

  ~Simulation();

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

  void init(std::size_t time_step);

 private:
  double _dx, _dy, _dz;
  double _cfl;
  int _num_thread;
  Divider::Type _divider_type;
  std::barrier<> _barrier;

  std::vector<std::shared_ptr<xfdtd::Object>> _objects;
  std::vector<std::shared_ptr<WaveformSource>> _waveform_sources;
  std::vector<std::shared_ptr<Boundary>> _boundaries;
  std::vector<std::shared_ptr<Monitor>> _monitors;
  std::vector<std::shared_ptr<Network>> _networks;
  std::vector<std::shared_ptr<NFFFT>> _nfffts;

  std::chrono::high_resolution_clock::time_point _start_time;
  std::chrono::high_resolution_clock::time_point _end_time;

  std::shared_ptr<GridSpace> _grid_space;
  std::shared_ptr<CalculationParam> _calculation_param;
  std::shared_ptr<EMF> _emf;

  std::vector<std::unique_ptr<Domain>> _domains;

  void generateDomain();

  void generateGridSpace();

  void correctMaterialSpace();

  void correctUpdateCoefficient();

  std::unique_ptr<Updator> makeUpdator(const Divider::IndexTask& task);

  void updateE();

  void updateH();

  void correctE();

  void correctH();

  void record();

  void printRunInfo();
};

}  // namespace xfdtd

#endif  // _XFDTD_LIB_SIMULATION_H_
