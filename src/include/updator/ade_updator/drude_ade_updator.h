#ifndef __XFDTD_CORE_UNIFORM_DISPERSIVE_UPDATOR_H__
#define __XFDTD_CORE_UNIFORM_DISPERSIVE_UPDATOR_H__

#include <xfdtd/material/ade_method/drude_ade_method.h>

#include <memory>

#include "updator/ade_updator/ade_updator.h"

namespace xfdtd {

class TemplateADEUpdateScheme;

class DrudeADEUpdator : public ADEUpdator {
  friend class TemplateADEUpdateScheme;

 public:
  DrudeADEUpdator(std::shared_ptr<const GridSpace> grid_space,
                  std::shared_ptr<CalculationParam> calculation_param,
                  std::shared_ptr<EMF> emf, IndexTask task,
                  std::shared_ptr<ADEMethodStorage> ade_method_storage);

 private:
  template <Axis::XYZ xyz>
  auto updateJ(Index i, Index j, Index k, const Real e_next,
               const Real e_cur) -> void;

  template <Axis::XYZ xyz>
  auto recordEPrevious(Real e, Index i, Index j, Index k) -> void {}

  template <Axis::XYZ xyz>
  auto calculateJSum(Index i, Index j, Index k) -> Real;

  template <Axis::XYZ xyz>
  auto ePrevious(Index i, Index j, Index k) const -> Real {
    return 0.0;
  }

  auto coeffEPrev(Index i, Index j, Index k) const -> Real { return 0.0; }
};

class DrudeADEUpdator3D : public DrudeADEUpdator {
 public:
  using DrudeADEUpdator::DrudeADEUpdator;

  auto updateE() -> void override;
};

}  // namespace xfdtd

#endif  // __XFDTD_CORE_UNIFORM_DISPERSIVE_UPDATOR_H__
