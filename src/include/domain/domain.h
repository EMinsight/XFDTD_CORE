#ifndef _XFDTD_CORE_DOMAIN_H_
#define _XFDTD_CORE_DOMAIN_H_


#include <xfdtd/grid_space/grid_space.h>
#include <xfdtd/monitor/monitor.h>
#include <xfdtd/nffft/nffft.h>
#include <xfdtd/waveform_source/waveform_source.h>

#include <barrier>
#include <memory>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "corrector/corrector.h"
#include "updator/updator.h"

namespace xfdtd {

class Domain {
 public:
  Domain(std::size_t id, IndexTask task,
         std::shared_ptr<GridSpace> grid_space,
         std::shared_ptr<CalculationParam> calculation_param,
         std::shared_ptr<EMF> emf, std::unique_ptr<Updator> updator,
         std::vector<std::shared_ptr<WaveformSource>> waveform_sources,
         std::vector<std::unique_ptr<Corrector>> correctors,
         std::vector<std::shared_ptr<Monitor>> monitors,
         std::vector<std::shared_ptr<NFFFT>> nfffts, std::barrier<>& barrier,
         bool master = false);

  Domain(std::size_t id, IndexTask task,
         std::shared_ptr<GridSpace> grid_space,
         std::shared_ptr<CalculationParam> calculation_param,
         std::shared_ptr<EMF> emf, std::unique_ptr<Updator> updator,
         std::barrier<>& barrier, bool master = false)
      : _id{id},
        _task{task},
        _grid_space{std::move(grid_space)},
        _calculation_param{std::move(calculation_param)},
        _emf{std::move(emf)},
        _updator{std::move(updator)},
        _barrier{barrier},
        _master{master} {}

  ~Domain() = default;

  auto id() const { return _id; }

  auto task() const { return _task; }

  bool isCalculationDone() const;

  void run();

  void updateE();

  void updateH();

  void correctE();

  void correctH();

  void record();

  void nextStep();

  void threadSynchronize();

  void processSynchronize();

  void synchronize();

  bool isMaster() const { return _master; }

  bool containXNEdge() const { return _task.xRange().start() == 0; }

  bool containYNEdge() const { return _task.yRange().start() == 0; }

  bool containZNEdge() const { return _task.zRange().start() == 0; }

  bool nodeContainXNBoundary() const {
    return _grid_space->globalBox().origin().i() == 0;
  }

  bool nodeContainYNBoundary() const {
    return _grid_space->globalBox().origin().j() == 0;
  }

  bool nodeContainZNBoundary() const {
    return _grid_space->globalBox().origin().k() == 0;
  }

  bool nodeContainXPBoundary() const {
    return _grid_space->globalBox().end().i() ==
           _grid_space->globalGridSpace()->sizeX();
  }

  bool nodeContainYPBoundary() const {
    return _grid_space->globalBox().end().j() ==
           _grid_space->globalGridSpace()->sizeY();
  }

  bool nodeContainZPBoundary() const {
    return _grid_space->globalBox().end().k() ==
           _grid_space->globalGridSpace()->sizeZ();
  }

  auto toString() const -> std::string;

 protected:
  void exchangeH();

 private:
  std::size_t _id;
  IndexTask _task;
  std::shared_ptr<GridSpace> _grid_space;
  std::shared_ptr<CalculationParam> _calculation_param;
  std::shared_ptr<EMF> _emf;
  std::unique_ptr<Updator> _updator;
  std::vector<std::shared_ptr<WaveformSource>> _waveform_sources{};
  std::vector<std::unique_ptr<Corrector>> _correctors{};
  std::vector<std::shared_ptr<Monitor>> _monitors{};
  std::vector<std::shared_ptr<NFFFT>> _nfffts{};
  std::barrier<>& _barrier;
  bool _master = false;
};

}  // namespace xfdtd

#endif  // _XFDTD_CORE_DOMAIN_H_
