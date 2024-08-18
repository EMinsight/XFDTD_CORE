#ifndef __XFDTD_CORE_M_LOR_ADE_UPDATOR_H__
#define __XFDTD_CORE_M_LOR_ADE_UPDATOR_H__

#include <xfdtd/material/ade_method/m_lor_ade_method.h>

#include "updator/basic_updator.h"

namespace xfdtd {

class MLorentzUpdator : public BasicUpdator3D {
 public:
  MLorentzUpdator(
      std::shared_ptr<const GridSpace> grid_space,
      std::shared_ptr<CalculationParam> calculation_param,
      std::shared_ptr<EMF> emf, IndexTask task,
      std::shared_ptr<MLorentzADEMethodStorage> m_lor_ade_method_storage);

  auto updateE() -> void override;

  auto& storage() const { return _m_lor_ade_method_storage; }

  auto& storage() { return _m_lor_ade_method_storage; }

 private:
  std::shared_ptr<MLorentzADEMethodStorage> _m_lor_ade_method_storage{};

  template <Axis::XYZ xyz>
  auto updateE() -> void;

  template <Axis::XYZ xyz>
  auto updateJ(Index i, Index j, Index k, const Real e_next,
               const Real e_cur) -> void;

  template <Axis::XYZ xyz>
  auto recordEPrevious(Real e, Index i, Index j, Index k) -> void;

  template <Axis::XYZ xyz>
  auto calculateJSum(Index i, Index j, Index k) -> Real;

  template <Axis::XYZ xyz>
  auto ePrevious(Index i, Index j, Index k) const -> Real;

  auto coeffEPrev(Index i, Index j, Index k) const -> Real;
};

}  // namespace xfdtd

#endif  // __XFDTD_CORE_M_LOR_ADE_UPDATOR_H__
