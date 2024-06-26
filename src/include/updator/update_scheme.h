#ifndef _XFDTD_CORE_UPDATE_SCHEME_H_
#define _XFDTD_CORE_UPDATE_SCHEME_H_

#include "xfdtd/calculation_param/calculation_param.h"
#include "xfdtd/coordinate_system/coordinate_system.h"
#include "xfdtd/electromagnetic_field/electromagnetic_field.h"
namespace xfdtd {

/**
 * FDTD Scheme for updating the electric field.
 * `a,b,c` means the axis.
 * `c` = `a` cross `b`.
 * `p,q` means the different direction.
 * `q` is the behind direction of `p`.
 *
 * @param cece The coefficient cece.
 * @param ec The current value of the electric field.
 * @param cecha The coefficient cecha.
 * @param ha_p The current value of the magnetic field ha_p.
 * @param ha_q The current value of the magnetic field ha_q.
 * @param cechb The coefficient cechb.
 * @param hb_p The current value of the magnetic field hb_p.
 * @param hb_q The current value of the magnetic field hb_q.
 * @return The next time step value of the electric field.
 */
template <typename T>
inline auto eNext(const T& cece, T ec, const T& cecha, T ha_p, T ha_q,
                  const T& cechb, T hb_p, T hb_q) {
  auto&& result = cece * ec + cecha * (ha_p - ha_q) + cechb * (hb_p - hb_q);
  return result;
}

/**
 * FDTD Scheme for updating the magnetic field.
 * `a,b,c` means the axis.
 * `c` = `a` cross `b`.
 * `p,q` means the different direction.
 * `q` is the behind direction of `p`.
 *
 * @param chch The coefficient chch.
 * @param hc The current value of the magnetic field.
 * @param chcea The coefficient chcea.
 * @param ea_p The current value of the electric field ea_p.
 * @param ea_q The current value of the electric field ea_q.
 * @param chceb The coefficient chceb.
 * @param eb_p The current value of the electric field eb_p.
 * @param eb_q The current value of the electric field eb_q.
 * @return The next time step value of the magnetic field.
 */
template <typename T>
inline auto hNextDeparted(const T& chch, T hc, const T& chcea, T ea_p, T ea_q,
                          const T& chceb, T eb_p, T eb_q) {
  auto&& result = chch * hc + chcea * (ea_p - ea_q) + chceb * (eb_p - eb_q);
  return result;
}

template <typename T>
inline auto hNext(const T& chch, T hc, const T& chcea, T ea_p, T ea_q,
                  const T& chceb, T eb_p, T eb_q) {
  auto&& result = chch * hc + chcea * (ea_q - ea_p) + chceb * (eb_q - eb_p);
  return result;
}

template <typename EMF::Attribute attribute, Axis::XYZ xyz, typename Size>
inline auto update(EMF& emf, FDTDUpdateCoefficient& update_coefficient,
                   const Size is, const Size ie, const Size js, const Size je,
                   const Size ks, const Size ke) {
  constexpr auto dual_attribute = EMF::dualAttribute(attribute);
  constexpr auto xzy_a = Axis::tangentialAAxis<xyz>();
  constexpr auto xzy_b = Axis::tangentialBAxis<xyz>();

  const auto& cfcf = update_coefficient.coeff<attribute, xyz>();
  const auto& cf_a =
      update_coefficient.coeff<attribute, xyz, dual_attribute, xzy_a>();
  const auto& cf_b =
      update_coefficient.coeff<attribute, xyz, dual_attribute, xzy_b>();

  auto&& field = emf.field<attribute, xyz>();
  const auto& field_a = emf.field<dual_attribute, xzy_a>();
  const auto& field_b = emf.field<dual_attribute, xzy_b>();

  constexpr Size offset = attribute == EMF::Attribute::E ? -1 : 1;
  for (Size i = is; i < ie; ++i) {
    for (Size j = js; j < je; ++j) {
      for (Size k = ks; k < ke; ++k) {
        auto [a, b, c] = transform::xYZToABC<Index, xyz>(i, j, k);
        auto b_1 = b + offset;
        auto a_1 = a + offset;

        auto [i_a, j_a, k_a] = transform::aBCToXYZ<Index, xyz>(a, b_1, c);

        auto [i_b, j_b, k_b] = transform::aBCToXYZ<Index, xyz>(a_1, b, c);

        if constexpr (attribute == EMF::Attribute::E) {
          field(i, j, k) =
              eNext(cfcf(i, j, k), field(i, j, k), cf_a(i, j, k),
                    field_a(i, j, k), field_a(i_a, j_a, k_a), cf_b(i, j, k),
                    field_b(i, j, k), field_b(i_b, j_b, k_b));
        } else {
          field(i, j, k) =
              hNext(cfcf(i, j, k), field(i, j, k), cf_a(i, j, k),
                    field_a(i, j, k), field_a(i_a, j_a, k_a), cf_b(i, j, k),
                    field_b(i, j, k), field_b(i_b, j_b, k_b));
        }
      }
    }
  }
}

}  // namespace xfdtd

#endif  // _XFDTD_CORE_UPDATE_SCHEME_H_
