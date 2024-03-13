#ifndef _XFDTD_CORE_LUMPED_ELEMENT_CORRECTOR_H_
#define _XFDTD_CORE_LUMPED_ELEMENT_CORRECTOR_H_

#include <xfdtd/divider/divider.h>

#include <memory>
#include <utility>

#include "corrector/corrector.h"
#include "xfdtd/calculation_param/calculation_param.h"

namespace xfdtd {

class LumpedElementCorrector : public Corrector {
 public:
  LumpedElementCorrector(Divider::IndexTask task, Divider::IndexTask local_task,
                         std::shared_ptr<CalculationParam> calculation_param,
                         xt::xarray<double>& e_field)
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
  Divider::IndexTask _task, _local_task;
  std::shared_ptr<CalculationParam> _calculation_param;
  xt::xarray<double>& _e_field;
};

class VoltageSourceCorrector : public LumpedElementCorrector {
 public:
  VoltageSourceCorrector(Divider::IndexTask task, Divider::IndexTask local_task,
                         std::shared_ptr<CalculationParam> calculation_param,
                         xt::xarray<double>& e_field,
                         const xt::xarray<double>& coeff_v,
                         const xt::xarray<double>& waveform)
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
  const xt::xarray<double>& _coeff_v;
  const xt::xarray<double>& _waveform;
};

class CurrentSourceCorrector : public LumpedElementCorrector {
 public:
  CurrentSourceCorrector(Divider::IndexTask task, Divider::IndexTask local_task,
                         std::shared_ptr<CalculationParam> calculation_param,
                         xt::xarray<double>& e_field,
                         const xt::xarray<double>& coeff_i,
                         const xt::xarray<double>& waveform)
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
  const xt::xarray<double>& _coeff_i;
  const xt::xarray<double>& _waveform;
};

class InductorCorrector : public LumpedElementCorrector {
 public:
 public:
  InductorCorrector(Divider::IndexTask task, Divider::IndexTask local_task,
                    std::shared_ptr<CalculationParam> calculation_param,
                    xt::xarray<double>& e_field, xt::xarray<double>& j,
                    const xt::xarray<double>& cecjc,
                    const xt::xarray<double>& cjcec)
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
  xt::xarray<double>& _j;
  const xt::xarray<double>&_cecjc, &_cjcec;
};

}  // namespace xfdtd

#endif  // _XFDTD_CORE_LUMPED_ELEMENT_CORRECTOR_H_
