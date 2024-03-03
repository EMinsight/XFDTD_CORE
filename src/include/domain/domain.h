#ifndef _XFDTD_LIB_DOMAIN_H_
#define _XFDTD_LIB_DOMAIN_H_

#include <barrier>
#include <memory>
#include <utility>
#include <vector>

#include "corrector/corrector.h"
#include "divider/divider.h"
#include "updator/updator.h"
#include "xfdtd/grid_space/grid_space.h"
#include "xfdtd/monitor/monitor.h"
#include "xfdtd/nffft/nffft.h"
#include "xfdtd/waveform_source/waveform_source.h"

namespace xfdtd {

class Domain {
 public:
  Domain(std::size_t id, Divider::IndexTask task,
         std::shared_ptr<GridSpace> grid_space,
         std::shared_ptr<CalculationParam> calculation_param,
         std::shared_ptr<EMF> emf, std::unique_ptr<Updator> updator,
         std::vector<std::shared_ptr<WaveformSource>> waveform_sources,
         std::vector<std::unique_ptr<Corrector>> correctors,
         std::vector<std::shared_ptr<Monitor>> monitors,
         std::vector<std::shared_ptr<NFFFT>> nfffts, std::barrier<>& barrier,
         bool master = false)
      : _id{id},
        _task{task},
        _grid_space{std::move(grid_space)},
        _calculation_param{std::move(calculation_param)},
        _emf{std::move(emf)},
        _updator{std::move(updator)},
        _waveform_sources{std::move(waveform_sources)},
        _correctors{std::move(correctors)},
        _monitors{std::move(monitors)},
        _nfffts{std::move(nfffts)},
        _barrier{barrier},
        _master{master} {
    std::cout << "Domain: " << _id << " " << _task.toString() << '\n';
    std::cout << _updator->toString() << '\n' << '\n';
  }

  ~Domain() = default;

  auto id() const { return _id; }

  bool isCalculationDone() const;

  void run();

  void updateE();

  void updateH();

  void correctE();

  void correctH();

  void synchronize();

  void communicate();

  void record();

  void nextStep();

  bool isMaster() const { return _master; }

  bool containXNEdge() const { return _task.xRange().start() == 0; }

  bool containYNEdge() const { return _task.yRange().start() == 0; }

  bool containZNEdge() const { return _task.zRange().start() == 0; }

 protected:
  void communicateForH();

  void communicateForE();

 private:
  std::size_t _id;
  Divider::IndexTask _task;
  std::shared_ptr<GridSpace> _grid_space;
  std::shared_ptr<CalculationParam> _calculation_param;
  std::shared_ptr<EMF> _emf;
  std::unique_ptr<Updator> _updator;
  std::vector<std::shared_ptr<WaveformSource>> _waveform_sources;
  std::vector<std::unique_ptr<Corrector>> _correctors;
  std::vector<std::shared_ptr<Monitor>> _monitors;
  std::vector<std::shared_ptr<NFFFT>> _nfffts;
  std::barrier<>& _barrier;
  bool _master = false;

  xt::xarray<double> makeHyBufferForEx();

  xt::xarray<double> makeHzBufferForEx();

  xt::xarray<double> makeHzBufferForEy();

  xt::xarray<double> makeHxBufferForEy();

  xt::xarray<double> makeHxBufferForEz();

  xt::xarray<double> makeHyBufferForEz();
};

}  // namespace xfdtd

#endif  // _XFDTD_LIB_DOMAIN_H_
