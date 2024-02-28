#ifndef _XFDTD_LIB_DOMAIN_H_
#define _XFDTD_LIB_DOMAIN_H_

#include <memory>
#include <vector>

#include "corrector/corrector.h"
#include "xfdtd/grid_space/grid_space.h"

namespace xfdtd {

class GridSpace;
class CalculationParam;
class EMF;
class Updator;
// class Corrector;
class Object;
class WaveformSource;

class Domain {
 public:
 private:
  GridBox _grid_box;
  std::shared_ptr<CalculationParam> _calculation_param;
  std::shared_ptr<EMF> _emf;
  std::unique_ptr<Updator> _updator;
  std::vector<Corrector> _corrector;
  std::vector<std::shared_ptr<Object>> _objects;
  std::vector<std::shared_ptr<WaveformSource>> _waveform_sources;

  bool isCalculationDone();

  void updateE();

  void updateH();

  void correctE();

  void correctH();

  void communicate();

  void record();
};

}  // namespace xfdtd

#endif  // _XFDTD_LIB_DOMAIN_H_
