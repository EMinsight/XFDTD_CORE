#ifndef _XFDTD_CORE_LUMPED_ELEMENT_CORRECTOR_H_
#define _XFDTD_CORE_LUMPED_ELEMENT_CORRECTOR_H_

#include <xfdtd/common/index_task.h>
#include <xfdtd/common/type_define.h>

#include <memory>
#include <utility>

#include "corrector/corrector.h"
#include "xfdtd/calculation_param/calculation_param.h"

namespace xfdtd {

class LumpedElementCorrector : public Corrector {
 public:
  LumpedElementCorrector(IndexTask task, IndexTask local_task,
                         std::shared_ptr<CalculationParam> calculation_param,
                         Array3D<Real>& e_field)
      : _task{task},
        _local_task{local_task},
        _calculation_param{std::move(calculation_param)},
        _e_field{e_field} {}

  ~LumpedElementCorrector() override = default;

  std::string toString() const override {
    std::stringstream ss;
    ss << "LumpedElementCorrector: ";
    ss << "task: " << _task.toString() << ", ";
    ss << "local_task: " << _local_task.toString() << ", ";
    return ss.str();
  }

 protected:
  IndexTask _task, _local_task;
  std::shared_ptr<CalculationParam> _calculation_param;
  Array3D<Real>& _e_field;
};

class VoltageSourceCorrector : public LumpedElementCorrector {
 public:
  VoltageSourceCorrector(IndexTask task, IndexTask local_task,
                         std::shared_ptr<CalculationParam> calculation_param,
                         Array3D<Real>& e_field, const Array3D<Real>& coeff_v,
                         const Array1D<Real>& waveform)
      : LumpedElementCorrector{task, local_task, std::move(calculation_param),
                               e_field},
        _coeff_v{coeff_v},
        _waveform{waveform} {}

  ~VoltageSourceCorrector() override = default;

  void correctE() override;

  void correctH() override;

  std::string toString() const override {
    std::stringstream ss;
    ss << "VoltageSourceCorrector: ";
    ss << LumpedElementCorrector::toString();
    return ss.str();
  }

 private:
  const Array3D<Real>& _coeff_v;
  const Array1D<Real>& _waveform;
};

class CurrentSourceCorrector : public LumpedElementCorrector {
 public:
  CurrentSourceCorrector(IndexTask task, IndexTask local_task,
                         std::shared_ptr<CalculationParam> calculation_param,
                         Array3D<Real>& e_field, const Array3D<Real>& coeff_i,
                         const Array1D<Real>& waveform)
      : LumpedElementCorrector{task, local_task, std::move(calculation_param),
                               e_field},
        _coeff_i{coeff_i},
        _waveform{waveform} {}

  ~CurrentSourceCorrector() override = default;

  void correctE() override;

  void correctH() override;

  std::string toString() const override {
    std::stringstream ss;
    ss << "CurrentSourceCorrector: ";
    ss << LumpedElementCorrector::toString() << "\n";
    return ss.str();
  }

 private:
  const Array3D<Real>& _coeff_i;
  const Array1D<Real>& _waveform;
};

class InductorCorrector : public LumpedElementCorrector {
 public:
 public:
  InductorCorrector(IndexTask task, IndexTask local_task,
                    std::shared_ptr<CalculationParam> calculation_param,
                    Array3D<Real>& e_field, Array3D<Real>& j,
                    const Array3D<Real>& cecjc, const Array3D<Real>& cjcec)
      : LumpedElementCorrector{task, local_task, std::move(calculation_param),
                               e_field},
        _j{j},
        _cecjc{cecjc},
        _cjcec{cjcec} {}

  void correctE() override;

  void correctH() override;

  std::string toString() const override {
    std::stringstream ss;
    ss << "InductorCorrector: ";
    ss << LumpedElementCorrector::toString() << "\n";
    return ss.str();
  }

 private:
  Array3D<Real>& _j;
  const Array3D<Real>&_cecjc, &_cjcec;
};

}  // namespace xfdtd

#endif  // _XFDTD_CORE_LUMPED_ELEMENT_CORRECTOR_H_
