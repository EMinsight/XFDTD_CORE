#ifndef __XFDTD_CORE_ADE_METHOD_H__
#define __XFDTD_CORE_ADE_METHOD_H__

#include <xfdtd/common/type_define.h>
#include <xfdtd/coordinate_system/coordinate_system.h>
#include <xfdtd/material/dispersive_material.h>

namespace xfdtd {

class ADEMethodStorage {
 public:
  ADEMethodStorage(Index num_pole);

  virtual ~ADEMethodStorage() = default;

  ADEMethodStorage(const ADEMethodStorage&) = default;

  ADEMethodStorage(ADEMethodStorage&&) noexcept = default;

  auto operator=(const ADEMethodStorage&) -> ADEMethodStorage& = default;

  auto operator=(ADEMethodStorage&&) noexcept -> ADEMethodStorage& = default;

  auto numPole() const { return _num_pole; }

  auto& coeffJJ() { return _coeff_j_j; }

  auto& coeffJJ() const { return _coeff_j_j; }

  auto& coeffJJPrev() { return _coeff_j_j_p; }

  auto& coeffJJPrev() const { return _coeff_j_j_p; }

  auto& coeffJENext() { return _coeff_j_e_n; }

  auto& coeffJENext() const { return _coeff_j_e_n; }

  auto& coeffJE() { return _coeff_j_e; }

  auto& coeffJE() const { return _coeff_j_e; }

  auto& coeffJEPrev() { return _coeff_j_e_p; }

  auto& coeffJEPrev() const { return _coeff_j_e_p; }

  auto& coeffJSumJ() { return _coeff_j_sum_j; }

  auto& coeffJSumJ() const { return _coeff_j_sum_j; }

  auto& coeffJSumJPrev() { return _coeff_j_sum_j_p; }

  auto& coeffJSumJPrev() const { return _coeff_j_sum_j_p; }

  auto& coeffEJSum() { return _coeff_e_j_sum; }

  auto& coeffEJSum() const { return _coeff_e_j_sum; }

  auto& coeffEEPrev() { return _coeff_e_e_p; }

  auto& coeffEEPrev() const { return _coeff_e_e_p; }

  auto& exPrev() { return _ex_prev; }

  auto& exPrev() const { return _ex_prev; }

  auto& eyPrev() { return _ey_prev; }

  auto& eyPrev() const { return _ey_prev; }

  auto& ezPrev() { return _ez_prev; }

  auto& ezPrev() const { return _ez_prev; }

  auto& jxArr() { return _jx_arr; }

  auto& jxArr() const { return _jx_arr; }

  auto& jyArr() { return _jy_arr; }

  auto& jyArr() const { return _jy_arr; }

  auto& jzArr() { return _jz_arr; }

  auto& jzArr() const { return _jz_arr; }

  auto& jxPrevArr() { return _jx_prev_arr; }

  auto& jxPrevArr() const { return _jx_prev_arr; }

  auto& jyPrevArr() { return _jy_prev_arr; }

  auto& jyPrevArr() const { return _jy_prev_arr; }

  auto& jzPrevArr() { return _jz_prev_arr; }

  auto& jzPrevArr() const { return _jz_prev_arr; }

  template <Axis::XYZ xyz>
  auto& ePrevious() {
    return const_cast<Array3D<Real>&>(
        static_cast<const ADEMethodStorage*>(this)->ePrevious<xyz>());
  }

  template <Axis::XYZ xyz>
  auto& ePrevious() const {
    if constexpr (xyz == Axis::XYZ::X) {
      return exPrev();
    } else if constexpr (xyz == Axis::XYZ::Y) {
      return eyPrev();
    } else if constexpr (xyz == Axis::XYZ::Z) {
      return ezPrev();
    }
  }

  template <Axis::XYZ xyz>
  auto& jArr() {
    return const_cast<Array4D<Real>&>(
        static_cast<const ADEMethodStorage*>(this)->jArr<xyz>());
  }

  template <Axis::XYZ xyz>
  auto& jArr() const {
    if constexpr (xyz == Axis::XYZ::X) {
      return jxArr();
    } else if constexpr (xyz == Axis::XYZ::Y) {
      return jyArr();
    } else if constexpr (xyz == Axis::XYZ::Z) {
      return jzArr();
    }
  }

  template <Axis::XYZ xyz>
  auto& jPrevArr() {
    return const_cast<Array4D<Real>&>(
        static_cast<const ADEMethodStorage*>(this)->jPrevArr<xyz>());
  }

  template <Axis::XYZ xyz>
  auto& jPrevArr() const {
    if constexpr (xyz == Axis::XYZ::X) {
      return jxPrevArr();
    } else if constexpr (xyz == Axis::XYZ::Y) {
      return jyPrevArr();
    } else if constexpr (xyz == Axis::XYZ::Z) {
      return jzPrevArr();
    }
  }

  virtual auto correctCoeff(
      Index i, Index j, Index k,
      const LinearDispersiveMaterial& linear_dispersive_material,
      const std::shared_ptr<const GridSpace>& grid_space,
      const std::shared_ptr<CalculationParam>& calculation_param) -> void = 0;

 protected:
  Index _num_pole{};
  Array4D<Real> _coeff_j_j{}, _coeff_j_j_p{}, _coeff_j_e_n{}, _coeff_j_e{},
      _coeff_j_e_p{}, _coeff_j_sum_j{}, _coeff_j_sum_j_p{};
  Array3D<Real> _coeff_e_j_sum{}, _coeff_e_e_p{};

  Array3D<Real> _ex_prev{}, _ey_prev{}, _ez_prev{};
  Array4D<Real> _jx_arr{}, _jy_arr{}, _jz_arr{}, _jx_prev_arr{}, _jy_prev_arr{},
      _jz_prev_arr{};
};

}  // namespace xfdtd

#endif  // __XFDTD_CORE_ADE_METHOD_H__
