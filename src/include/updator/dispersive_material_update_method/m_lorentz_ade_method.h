#ifndef __XFDTD_CORE_M_LORENTZ_J_E_METHOD_H__
#define __XFDTD_CORE_M_LORENTZ_J_E_METHOD_H__

#include <xfdtd/common/type_define.h>
#include <xfdtd/material/dispersive_material.h>

#include "updator/dispersive_material_update_method/dispersive_material_update_method.h"

namespace xfdtd {

class MLorentzADEMethod : public LinearDispersiveMaterialUpdateMethod {
 public:
  static auto a(const MLorentzEqDecision& m_lorentz_equation, Real dt)
      -> Array1D<Real>;

  static auto b(const MLorentzEqDecision& m_lorentz_equation, Real dt)
      -> Array1D<Real>;

  static auto c(const MLorentzEqDecision& m_lorentz_equation, Real dt)
      -> Array1D<Real>;

  static auto d(const MLorentzEqDecision& m_lorentz_equation, Real dt)
      -> Array1D<Real>;

  static auto e(const MLorentzEqDecision& m_lorentz_equation, Real dt)
      -> Array1D<Real>;

  static auto f(const MLorentzEqDecision& m_lorentz_equation, Real dt)
      -> Array1D<Real>;

 public:
  MLorentzADEMethod(Real epsilon_inf,
                    std::shared_ptr<MLorentzEqDecision> m_lorentz_equation);

  ~MLorentzADEMethod() override = default;

  auto clone() const
      -> std::unique_ptr<LinearDispersiveMaterialUpdateMethod> override {
    return std::make_unique<MLorentzADEMethod>(*this);
  }

  auto init(Real dt) -> void override;

  auto correctCoeff(Index i, Index j, Index k, const GridSpace* grid_space,
                    CalculationParam* calculation_param) const -> void override;

  auto initUpdate(const GridSpace* grid_space,
                  std::shared_ptr<const CalculationParam> calculation_param,
                  std::shared_ptr<EMF> emf, Index m_index,
                  const IndexTask& task) -> void override;

  auto updateEx(Index i, Index j, Index k) -> void override;

  auto updateEy(Index i, Index j, Index k) -> void override;

  auto updateEz(Index i, Index j, Index k) -> void override;

 private:
  std::shared_ptr<MLorentzEqDecision> _m_lorentz_equation;

  Array1D<Real> _a, _b, _c, _d, _e, _f;

  Array1D<Real> _coeff_j_e_n, _coeff_j_e_c, _coeff_j_e_p, _coeff_j_j_c,
      _coeff_j_j_p;

  Array1D<Real> _coeff_j_sum_j_c, _coeff_j_sum_j_p;

  Real _coeff_e_e_p, _coeff_e_j_sum;

  Array3D<Real> _ex_prev, _ey_prev, _ez_prev;

  Array4D<Real> _jx, _jy, _jz, _jx_prev, _jy_prev, _jz_prev;

  auto calculateJSum(Index i, Index j, Index k, const Array4D<Real>& j_arr,
                     const Array4D<Real>& j_prev_arr) -> Real;

  auto updateJ(Index i, Index j, Index k, Real e_next, Real e_cur, Real e_prev,
               Array4D<Real>& j_arr, Array4D<Real>& j_prev_arr) -> void;

  auto checkStabilityCondition(Real dt) -> void;
};

}  // namespace xfdtd

#endif  // __XFDTD_CORE_M_LORENTZ_J_E_METHOD_H__
