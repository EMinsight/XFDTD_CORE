#include "updator/ade_updator/drude_ade_updator.h"

#include "updator/ade_updator/template_ade_update_scheme.h"

namespace xfdtd {

DrudeADEUpdator::DrudeADEUpdator(
    std::shared_ptr<const GridSpace> grid_space,
    std::shared_ptr<CalculationParam> calculation_param,
    std::shared_ptr<EMF> emf, IndexTask task,
    std::shared_ptr<ADEMethodStorage> ade_method_storage)
    : ADEUpdator{grid_space, calculation_param, emf, task, ade_method_storage} {
}

auto DrudeADEUpdator3D::updateE() -> void {
  TemplateADEUpdateScheme::updateE<DrudeADEUpdator, Axis::XYZ::X>(this);
  TemplateADEUpdateScheme::updateE<DrudeADEUpdator, Axis::XYZ::Y>(this);
  TemplateADEUpdateScheme::updateE<DrudeADEUpdator, Axis::XYZ::Z>(this);
}

auto drudeCalculateJSum(Index i, Index j, Index k, Index num_p,
                        const Array4D<Real>& j_arr,
                        const Array4D<Real>& coeff_j_sum_j) -> Real {
  Real sum{0};
  for (Index p{0}; p < num_p; ++p) {
    sum += coeff_j_sum_j(i, j, k, p) * j_arr(i, j, k, p);
  }
  return sum;
}

auto drudeUpdateJ(Index i, Index j, Index k, Index num_p, Array4D<Real>& j_arr,
                  const Array4D<Real>& coeff_j_j,
                  const Array4D<Real>& coeff_j_e, Real e_next,
                  Real e_cur) -> void {
  for (Index p{0}; p < num_p; ++p) {
    j_arr(i, j, k, p) = coeff_j_j(i, j, k, p) * j_arr(i, j, k, p) +
                        coeff_j_e(i, j, k, p) * (e_next + e_cur);
  }
}

template <Axis::XYZ xyz>
auto DrudeADEUpdator::calculateJSum(Index i, Index j, Index k) -> Real {
  const auto& drude_ade_method_storage = this->storage();
  auto num_p = this->storage()->numPole();
  const auto& j_arr = this->storage()->jArr<xyz>();
  const auto& coeff_j_sum_j = this->storage()->coeffJSumJ();
  return drudeCalculateJSum(i, j, k, num_p, j_arr, coeff_j_sum_j);
}

template <Axis::XYZ xyz>
auto DrudeADEUpdator::updateJ(Index i, Index j, Index k, const Real e_next,
                              const Real e_cur) -> void {
  auto num_p = this->storage()->numPole();
  auto& j_arr = this->storage()->jArr<xyz>();
  const auto& coeff_j_j = this->storage()->coeffJJ();
  const auto& coeff_j_e = this->storage()->coeffJE();
  drudeUpdateJ(i, j, k, num_p, j_arr, coeff_j_j, coeff_j_e, e_next, e_cur);
  return;
}

// explicit instantiation
template auto DrudeADEUpdator::calculateJSum<Axis::XYZ::X>(Index i, Index j,
                                                           Index k) -> Real;

template auto DrudeADEUpdator::calculateJSum<Axis::XYZ::Y>(Index i, Index j,
                                                           Index k) -> Real;

template auto DrudeADEUpdator::calculateJSum<Axis::XYZ::Z>(Index i, Index j,
                                                           Index k) -> Real;

template auto DrudeADEUpdator::updateJ<Axis::XYZ::X>(Index i, Index j, Index k,
                                                     const Real e_next,
                                                     const Real e_cur) -> void;

template auto DrudeADEUpdator::updateJ<Axis::XYZ::Y>(Index i, Index j, Index k,
                                                     const Real e_next,
                                                     const Real e_cur) -> void;

template auto DrudeADEUpdator::updateJ<Axis::XYZ::Z>(Index i, Index j, Index k,
                                                     const Real e_next,
                                                     const Real e_cur) -> void;

}  // namespace xfdtd
