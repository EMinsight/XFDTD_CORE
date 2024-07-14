#ifndef _XFDTD_CORE_WAVEFORM_SOURCE_CORRECTOR_H_
#define _XFDTD_CORE_WAVEFORM_SOURCE_CORRECTOR_H_

#include <xfdtd/calculation_param/calculation_param.h>
#include <xfdtd/common/index_task.h>
#include <xfdtd/common/type_define.h>
#include <xfdtd/electromagnetic_field/electromagnetic_field.h>
#include <xfdtd/grid_space/grid_space.h>

#include "corrector/corrector.h"

namespace xfdtd {

class WaveformSourceCorrector : public Corrector {
 public:
  WaveformSourceCorrector(IndexTask task, IndexTask local_task,
                          GridSpace* grid_space,
                          CalculationParam* calculation_param, EMF* emf,
                          const Array1D<Real>* waveform)
      : _task{task},
        _local_task{local_task},
        _grid_space{(grid_space)},
        _calculation_param{(calculation_param)},
        _emf{(emf)},
        _waveform{waveform} {}

  ~WaveformSourceCorrector() override = default;

  auto task() -> IndexTask { return _task; }

  auto localTask() -> IndexTask { return _local_task; }

  auto gridSpace() const { return _grid_space; }

  auto calculationParam() const { return _calculation_param; }

  auto emf() const { return _emf; }

  auto gridSpace() { return _grid_space; }

  auto calculationParam() { return _calculation_param; }

  auto emf() { return _emf; }

  auto waveform() -> const Array1D<Real>& { return *_waveform; }

  std::string toString() const override { return "WaveformSourceCorrector"; }

 private:
  IndexTask _task, _local_task;
  GridSpace* _grid_space;
  CalculationParam* _calculation_param;
  EMF* _emf;
  const Array1D<Real>* _waveform;
};

}  // namespace xfdtd

#endif  // _XFDTD_CORE_WAVEFORM_SOURCE_CORRECTOR_H_
