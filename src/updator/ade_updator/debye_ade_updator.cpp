#include "updator/ade_updator/debye_ade_updator.h"

#include <xfdtd/material/ade_method/ade_method.h>

#include "updator/ade_updator/ade_updator.h"
#include "updator/ade_updator/template_ade_update_scheme.h"
#include "updator/updator.h"

namespace xfdtd {

DebyeADEUpdator::DebyeADEUpdator(
    std::shared_ptr<const GridSpace> grid_space,
    std::shared_ptr<CalculationParam> calculation_param,
    std::shared_ptr<EMF> emf, IndexTask task,
    std::shared_ptr<ADEMethodStorage> storage)
    : ADEUpdator{grid_space, calculation_param, emf, task, storage} {}

auto DebyeADEUpdator3D::updateE() -> void {
  TemplateADEUpdateScheme::updateE<DebyeADEUpdator, Axis::XYZ::X>(this);
  TemplateADEUpdateScheme::updateE<DebyeADEUpdator, Axis::XYZ::Y>(this);
  TemplateADEUpdateScheme::updateE<DebyeADEUpdator, Axis::XYZ::Z>(this);
}

auto DebyeADEUpdator2D::updateE() -> void {
  TemplateADEUpdateScheme::updateE<DebyeADEUpdator, Axis::XYZ::Z>(this);
}

auto DebyeADEUpdator1D::updateE() -> void {
  // TemplateADEUpdateScheme::updateE<DebyeADEUpdator, Axis::XYZ::X>(this);
  throw XFDTDUpdatorException("1D DebyeADEUpdator is not implemented yet.");
}

static auto debyeCalculateJSum(Index i, Index j, Index k, Index num_p,
                               const Array4D<Real>& j_arr,
                               const Array4D<Real>& coeff_j_sum_j) -> Real {
  Real sum{0};
  for (Index p{0}; p < num_p; ++p) {
    sum += coeff_j_sum_j(i, j, k, p) * j_arr(i, j, k, p);
  }
  return sum;
}

static auto debyeUpdateJ(Index i, Index j, Index k, Index num_p,
                         Array4D<Real>& j_arr, const Array4D<Real>& coeff_j_j,
                         const Array4D<Real>& coeff_j_e, Real e_next,
                         Real e_cur) -> void {
  for (Index p{0}; p < num_p; ++p) {
    j_arr(i, j, k, p) = coeff_j_j(i, j, k, p) * j_arr(i, j, k, p) +
                        coeff_j_e(i, j, k, p) * (e_next - e_cur);
  }
}

template <Axis::XYZ xyz>
auto DebyeADEUpdator::updateJ(Index i, Index j, Index k, const Real e_next,
                              const Real e_cur) -> void {
  auto& j_arr = this->storage()->jArr<xyz>();
  const auto& coeff_j_j = this->storage()->coeffJJ();
  const auto& coeff_j_e = this->storage()->coeffJE();
  const auto num_p = this->storage()->numPole();
  debyeUpdateJ(i, j, k, num_p, j_arr, coeff_j_j, coeff_j_e, e_next, e_cur);
}

template <Axis::XYZ xyz>
auto DebyeADEUpdator::calculateJSum(Index i, Index j, Index k) -> Real {
  const auto& j_arr = this->storage()->jArr<xyz>();
  const auto& coeff_j_sum_j = this->storage()->coeffJSumJ();
  const auto num_p = this->storage()->numPole();
  return debyeCalculateJSum(i, j, k, num_p, j_arr, coeff_j_sum_j);
}

// explicit instantiation
template auto DebyeADEUpdator::calculateJSum<Axis::XYZ::X>(Index i, Index j,
                                                           Index k) -> Real;
template auto DebyeADEUpdator::calculateJSum<Axis::XYZ::Y>(Index i, Index j,
                                                           Index k) -> Real;
template auto DebyeADEUpdator::calculateJSum<Axis::XYZ::Z>(Index i, Index j,
                                                           Index k) -> Real;

template auto DebyeADEUpdator::updateJ<Axis::XYZ::X>(Index i, Index j, Index k,
                                                     const Real e_next,
                                                     const Real e_cur) -> void;
template auto DebyeADEUpdator::updateJ<Axis::XYZ::Y>(Index i, Index j, Index k,
                                                     const Real e_next,
                                                     const Real e_cur) -> void;
template auto DebyeADEUpdator::updateJ<Axis::XYZ::Z>(Index i, Index j, Index k,
                                                     const Real e_next,
                                                     const Real e_cur) -> void;

}  // namespace xfdtd
