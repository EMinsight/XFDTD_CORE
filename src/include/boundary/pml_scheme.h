#ifndef _XFDTD_CORE_PML_SCHEME_H_
#define _XFDTD_CORE_PML_SCHEME_H_

namespace xfdtd {

template <typename T>
inline auto correctPML(T& field, T& psi, const T& coeff_a, const T& coeff_b,
                       const T& field_p, const T& field_q, const T& c_psi) {
  psi = coeff_b * psi + coeff_a * (field_p - field_q);
  field += c_psi * psi;
}

}  // namespace xfdtd

#endif  // _XFDTD_CORE_PML_SCHEME_H_
