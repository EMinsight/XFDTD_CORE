#ifndef __XFDTD_CORE_DEBYE_ADE_UPDATOR_H__
#define __XFDTD_CORE_DEBYE_ADE_UPDATOR_H__

#include <xfdtd/material/ade_method/debye_ade_method.h>

#include "updator/basic_updator.h"

namespace xfdtd {

class DebyeADEUpdator : public BasicUpdator3D {
 public:
  DebyeADEUpdator(
      std::shared_ptr<const GridSpace> grid_space,
      std::shared_ptr<CalculationParam> calculation_param,
      std::shared_ptr<EMF> emf, IndexTask task,
      std::shared_ptr<DebyeADEMethodStorage> debye_ade_method_storage);

  auto updateE() -> void override;

  auto& storage() const { return _storage; }

  auto& storage() { return _storage; }

 private:
  std::shared_ptr<DebyeADEMethodStorage> _storage{};

  template <Axis::XYZ xyz>
  auto updateE() -> void;

  auto recordPrevE(Array3D<Real>& e_prev, const Array3D<Real>& e, Index i,
                   Index j, Index k) -> void;

  auto recordPrevJ(Array3D<Real>& j_prev, const Array3D<Real>& e, Index i,
                   Index j, Index k) -> void;

  template <Axis::XYZ xyz>
  auto calculateJSum(Index i, Index j, Index k) -> Real;

  template <Axis::XYZ xyz>
  auto updateJ(Index i, Index j, Index k, const Real e_next,
               const Real e_cur) -> void;
};

}  // namespace xfdtd

#endif  // __XFDTD_CORE_DEBYE_ADE_UPDATOR_H__
