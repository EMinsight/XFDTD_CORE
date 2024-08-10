#include "updator/ade_updator/m_lor_ade_updator.h"

#include "updator/update_scheme.h"

namespace xfdtd {

MLorentzUpdator::MLorentzUpdator(
    std::shared_ptr<const GridSpace> grid_space,
    std::shared_ptr<CalculationParam> calculation_param,
    std::shared_ptr<EMF> emf, IndexTask task,
    std::shared_ptr<MLorentzADEMethodStorage> m_lor_ade_method_storage)
    : BasicUpdator3D(grid_space, calculation_param, emf, task),
      _m_lor_ade_method_storage{m_lor_ade_method_storage} {}

auto MLorentzUpdator::updateE() -> void {
  this->updateE<Axis::XYZ::X>();
  this->updateE<Axis::XYZ::Y>();
  this->updateE<Axis::XYZ::Z>();
}

template <Axis::XYZ xyz>
auto MLorentzUpdator::updateE() -> void {
  const auto task = this->task();
  const auto& update_coefficient = this->calculationParam()->fdtdCoefficient();
  auto&& emf = this->emf();

  constexpr auto attribute = EMF::Attribute::E;
  constexpr auto dual_attribute = EMF::Attribute::H;
  constexpr auto xzy_a = Axis::tangentialAAxis<xyz>();  // !
  constexpr auto xzy_b = Axis::tangentialBAxis<xyz>();

  const auto& cfcf = update_coefficient->coeff<attribute, xyz>();
  const auto& cf_a =
      update_coefficient->coeff<attribute, xyz, dual_attribute, xzy_a>();
  const auto& cf_b =
      update_coefficient->coeff<attribute, xyz, dual_attribute, xzy_b>();

  auto&& field = emf->field<attribute, xyz>();
  const auto& field_a = emf->field<dual_attribute, xzy_a>();
  const auto& field_b = emf->field<dual_attribute, xzy_b>();

  auto is = task.xRange().start();
  auto ie = task.xRange().end();
  auto js = task.yRange().start();
  auto je = task.yRange().end();
  auto ks = task.zRange().start();
  auto ke = task.zRange().end();

  {
    auto [as, bs, cs] = transform::xYZToABC<decltype(is), xyz>(is, js, ks);
    as = (as == 0) ? 1 : as;
    bs = (bs == 0) ? 1 : bs;
    auto [i, j, k] = transform::aBCToXYZ<decltype(as), xyz>(as, bs, cs);
    is = i;
    js = j;
    ks = k;
  }

  constexpr Index offset = attribute == EMF::Attribute::E ? -1 : 1;

  auto&& storage = this->storage();

  const auto& coeff_e_j_sum = storage->coeffEJSum();

  for (Index i{is}; i < ie; ++i) {
    for (Index j{js}; j < je; ++j) {
      for (Index k{ks}; k < ke; ++k) {
        auto [a, b, c] = transform::xYZToABC<Index, xyz>(i, j, k);
        auto b_1 = b + offset;
        auto a_1 = a + offset;
        auto [i_a, j_a, k_a] = transform::aBCToXYZ<Index, xyz>(a, b_1, c);
        auto [i_b, j_b, k_b] = transform::aBCToXYZ<Index, xyz>(a_1, b, c);

        const auto j_sum = calculateJSum<xyz>(i, j, k);
        const auto e_prev = ePrevious<xyz>(i, j, k);
        const auto coeff_e_e_p = coeffEPrev(i, j, k);

        const auto e_cur = field(i, j, k);
        field(i, j, k) =
            coeff_e_e_p * e_prev +
            eNext(cfcf(i, j, k), field(i, j, k), cf_a(i, j, k),
                  field_a(i, j, k), field_a(i_a, j_a, k_a), cf_b(i, j, k),
                  field_b(i, j, k), field_b(i_b, j_b, k_b)) +
            coeff_e_j_sum(i, j, k) * j_sum;
        updateJ<xyz>(i, j, k, field(i, j, k), e_cur);

        recordEPrevious<xyz>(e_cur, i, j, k);
      }
    }
  }
}

template <Axis::XYZ xyz>
auto MLorentzUpdator::updateJ(Index i, Index j, Index k, const Real e_next,
                              const Real e_cur) -> void {
  const auto num_p = this->storage()->numPole();

  const auto& coeff_j_e_n = this->storage()->coeffJENext();
  const auto& coeff_j_e = this->storage()->coeffJE();
  const auto& coeff_j_e_p = this->storage()->coeffJEPrev();
  const auto& coeff_j_j = this->storage()->coeffJJ();
  const auto& coeff_j_j_p = this->storage()->coeffJJPrev();

  const auto e_prev = this->storage()->ePrevious<xyz>()(i, j, k);
  auto& j_arr = this->storage()->jArr<xyz>();
  auto& j_prev_arr = this->storage()->jPrevArr<xyz>();

  for (Index p{0}; p < num_p; ++p) {
    auto j_next = coeff_j_e_n(i, j, k, p) * e_next +
                  coeff_j_e(i, j, k, p) * e_cur +
                  coeff_j_e_p(i, j, k, p) * e_prev +
                  coeff_j_j(i, j, k, p) * j_arr(i, j, k, p) +
                  coeff_j_j_p(i, j, k, p) * j_prev_arr(i, j, k, p);
    j_prev_arr(i, j, k, p) = j_arr(i, j, k, p);
    j_arr(i, j, k, p) = j_next;
  }
}

template <Axis::XYZ xyz>
auto MLorentzUpdator::recordEPrevious(Real e, Index i, Index j,
                                      Index k) -> void {
  this->storage()->ePrevious<xyz>()(i, j, k) = e;
}

template <Axis::XYZ xyz>
auto MLorentzUpdator::calculateJSum(Index i, Index j, Index k) -> Real {
  Real sum = 0;
  const auto num_p = this->storage()->numPole();
  const auto& coeff_j_sum_j = this->storage()->coeffJSumJ();
  const auto& coeff_j_sum_j_p = this->storage()->coeffJSumJPrev();
  auto& j_arr = this->storage()->jArr<xyz>();
  auto& j_prev_arr = this->storage()->jPrevArr<xyz>();

  for (Index p{0}; p < this->storage()->numPole(); ++p) {
    sum += coeff_j_sum_j(i, j, k, p) * j_arr(i, j, k, p) +
           coeff_j_sum_j_p(i, j, k, p) * j_prev_arr(i, j, k, p);
  }

  return sum;
}

template <Axis::XYZ xyz>
auto MLorentzUpdator::ePrevious(Index i, Index j, Index k) const -> Real {
  return this->storage()->ePrevious<xyz>()(i, j, k);
}

auto MLorentzUpdator::coeffEPrev(Index i, Index j, Index k) const -> Real {
  return this->storage()->coeffEEPrev()(i, j, k);
}

// explicit instantiation
template auto MLorentzUpdator::updateE<Axis::XYZ::X>() -> void;
template auto MLorentzUpdator::updateE<Axis::XYZ::Y>() -> void;
template auto MLorentzUpdator::updateE<Axis::XYZ::Z>() -> void;

template auto MLorentzUpdator::calculateJSum<Axis::XYZ::X>(Index i, Index j,
                                                           Index k) -> Real;
template auto MLorentzUpdator::calculateJSum<Axis::XYZ::Y>(Index i, Index j,
                                                           Index k) -> Real;
template auto MLorentzUpdator::calculateJSum<Axis::XYZ::Z>(Index i, Index j,
                                                           Index k) -> Real;

template auto MLorentzUpdator::updateJ<Axis::XYZ::X>(Index i, Index j, Index k,
                                                     const Real e_next,
                                                     const Real e_cur) -> void;
template auto MLorentzUpdator::updateJ<Axis::XYZ::Y>(Index i, Index j, Index k,
                                                     const Real e_next,
                                                     const Real e_cur) -> void;
template auto MLorentzUpdator::updateJ<Axis::XYZ::Z>(Index i, Index j, Index k,
                                                     const Real e_next,
                                                     const Real e_cur) -> void;

template auto MLorentzUpdator::recordEPrevious<Axis::XYZ::X>(Real e, Index i,
                                                             Index j,
                                                             Index k) -> void;
template auto MLorentzUpdator::recordEPrevious<Axis::XYZ::Y>(Real e, Index i,
                                                             Index j,
                                                             Index k) -> void;
template auto MLorentzUpdator::recordEPrevious<Axis::XYZ::Z>(Real e, Index i,
                                                             Index j,
                                                             Index k) -> void;

template auto MLorentzUpdator::ePrevious<Axis::XYZ::X>(Index i, Index j,
                                                       Index k) const -> Real;
template auto MLorentzUpdator::ePrevious<Axis::XYZ::Y>(Index i, Index j,
                                                       Index k) const -> Real;
template auto MLorentzUpdator::ePrevious<Axis::XYZ::Z>(Index i, Index j,
                                                       Index k) const -> Real;

}  // namespace xfdtd
