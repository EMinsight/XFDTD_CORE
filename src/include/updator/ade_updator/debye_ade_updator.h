#ifndef __XFDTD_CORE_DEBYE_ADE_UPDATOR_H__
#define __XFDTD_CORE_DEBYE_ADE_UPDATOR_H__

#include <xfdtd/material/ade_method/ade_method.h>
#include <xfdtd/material/ade_method/debye_ade_method.h>

#include "updator/ade_updator/ade_updator.h"

namespace xfdtd {

class TemplateADEUpdateScheme;

class DebyeADEUpdator : public ADEUpdator {
  friend TemplateADEUpdateScheme;

 public:
  DebyeADEUpdator(std::shared_ptr<const GridSpace> grid_space,
                  std::shared_ptr<CalculationParam> calculation_param,
                  std::shared_ptr<EMF> emf, IndexTask task,
                  std::shared_ptr<ADEMethodStorage> debye_ade_method_storage);

 protected:
  template <Axis::XYZ xyz>
  auto updateE() -> void;

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

class DebyeADEUpdator3D : public DebyeADEUpdator {
 public:
  using DebyeADEUpdator::DebyeADEUpdator;

  auto updateE() -> void override;
};

class DebyeADEUpdator2D : public DebyeADEUpdator {
 public:
  using DebyeADEUpdator::DebyeADEUpdator;

  auto updateE() -> void override;
};

class DebyeADEUpdator1D : public DebyeADEUpdator {
 public:
  using DebyeADEUpdator::DebyeADEUpdator;

  auto updateE() -> void override;
};

}  // namespace xfdtd

#endif  // __XFDTD_CORE_DEBYE_ADE_UPDATOR_H__
