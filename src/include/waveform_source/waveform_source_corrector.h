#ifndef _XFDTD_CORE_WAVEFORM_SOURCE_CORRECTOR_H_
#define _XFDTD_CORE_WAVEFORM_SOURCE_CORRECTOR_H_

#include <xfdtd/calculation_param/calculation_param.h>
#include <xfdtd/common/index_task.h>
#include <xfdtd/common/type_define.h>
#include <xfdtd/electromagnetic_field/electromagnetic_field.h>
#include <xfdtd/grid_space/grid_space.h>

#include <memory>
#include <utility>

#include "corrector/corrector.h"

namespace xfdtd {

class WaveformSourceCorrector : public Corrector {
 public:
  WaveformSourceCorrector(IndexTask task, IndexTask local_task,
                          std::shared_ptr<GridSpace> grid_space,
                          std::shared_ptr<CalculationParam> calculation_param,
                          std::shared_ptr<EMF> emf,
                          const Array1D<Real>& waveform)
      : _task{task},
        _local_task{local_task},
        _grid_space{std::move(grid_space)},
        _calculation_param{std::move(calculation_param)},
        _emf{std::move(emf)},
        _waveform{waveform} {}

  ~WaveformSourceCorrector() override = default;

  auto task() -> IndexTask { return _task; }

  auto localTask() -> IndexTask { return _local_task; }

  auto gridSpace() -> std::shared_ptr<GridSpace> { return _grid_space; }

  auto calculationParam() -> std::shared_ptr<CalculationParam> {
    return _calculation_param;
  }

  auto emf() -> std::shared_ptr<EMF> { return _emf; }

  auto waveform() -> const Array1D<Real>& { return _waveform; }

  std::string toString() const override { return "WaveformSourceCorrector"; }

 protected:
  auto calculationParamPtr() { return _calculation_param.get(); }

  auto emfPtr() { return _emf.get(); }

 private:
  IndexTask _task, _local_task;
  std::shared_ptr<GridSpace> _grid_space;
  std::shared_ptr<CalculationParam> _calculation_param;
  std::shared_ptr<EMF> _emf;
  const Array1D<Real>& _waveform;
};

}  // namespace xfdtd

#endif  // _XFDTD_CORE_WAVEFORM_SOURCE_CORRECTOR_H_
