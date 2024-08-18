#ifndef __XFDTD_CORE_TEMPLATE_ADE_UPDATE_SCHEME_H__
#define __XFDTD_CORE_TEMPLATE_ADE_UPDATE_SCHEME_H__

#include <xfdtd/coordinate_system/coordinate_system.h>
#include <xfdtd/electromagnetic_field/electromagnetic_field.h>
#include <xfdtd/util/transform/abc_xyz.h>

#include "updator/update_scheme.h"

namespace xfdtd {

class TemplateADEUpdateScheme {
 public:
  template <typename ADEUpdator, Axis::XYZ xyz>
  static inline auto updateE(ADEUpdator* ade_updator) -> void {
    const auto task = ade_updator->task();
    const auto& update_coefficient =
        ade_updator->calculationParam()->fdtdCoefficient();
    auto&& emf = ade_updator->emf();

    constexpr auto attribute = EMF::Attribute::E;
    constexpr auto dual_attribute = EMF::Attribute::H;
    constexpr auto xzy_a = Axis::tangentialAAxis<xyz>();  // !
    constexpr auto xzy_b = Axis::tangentialBAxis<xyz>();

    const auto& cfcf = update_coefficient->template coeff<attribute, xyz>();
    const auto& cf_a =
        update_coefficient
            ->template coeff<attribute, xyz, dual_attribute, xzy_a>();
    const auto& cf_b =
        update_coefficient
            ->template coeff<attribute, xyz, dual_attribute, xzy_b>();

    auto&& field = emf->template field<attribute, xyz>();
    const auto& field_a = emf->template field<dual_attribute, xzy_a>();
    const auto& field_b = emf->template field<dual_attribute, xzy_b>();

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

    auto&& storage = ade_updator->storage();

    const auto& coeff_e_j_sum = storage->coeffEJSum();
    auto& field_prev = storage->template ePrevious<xyz>();

    for (Index i{is}; i < ie; ++i) {
      for (Index j{js}; j < je; ++j) {
        for (Index k{ks}; k < ke; ++k) {
          auto [a, b, c] = transform::xYZToABC<Index, xyz>(i, j, k);
          auto b_1 = b + offset;
          auto a_1 = a + offset;
          auto [i_a, j_a, k_a] = transform::aBCToXYZ<Index, xyz>(a, b_1, c);
          auto [i_b, j_b, k_b] = transform::aBCToXYZ<Index, xyz>(a_1, b, c);

          const auto j_sum = ade_updator->template calculateJSum<xyz>(i, j, k);
          const auto e_prev = ade_updator->template ePrevious<xyz>(i, j, k);
          const auto coeff_e_e_p = ade_updator->coeffEPrev(i, j, k);

          const auto e_cur = field(i, j, k);
          field(i, j, k) =
              coeff_e_e_p * e_prev +
              eNext(cfcf(i, j, k), field(i, j, k), cf_a(i, j, k),
                    field_a(i, j, k), field_a(i_a, j_a, k_a), cf_b(i, j, k),
                    field_b(i, j, k), field_b(i_b, j_b, k_b)) +
              coeff_e_j_sum.at(i, j, k) * j_sum;
          ade_updator->template updateJ<xyz>(i, j, k, field(i, j, k), e_cur);

          ade_updator->template recordEPrevious<xyz>(e_cur, i, j, k);
        }
      }
    }
  }
};

}  // namespace xfdtd

#endif  // __XFDTD_CORE_TEMPLATE_ADE_UPDATE_SCHEME_H__
