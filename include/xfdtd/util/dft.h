#ifndef _XFDTD_CORE_DFT_H_
#define _XFDTD_CORE_DFT_H_

#include <xfdtd/common/constant.h>
#include <xfdtd/common/type_define.h>

namespace xfdtd {

template <typename Arr>
inline Array1D<std::complex<Real>> dft(const Arr &time_domain_data, Real dt,
                                       const Arr &frequencies,
                                       Real time_shift = 0) {
  Array1D<std::complex<Real>> res =
      xt::empty<std::complex<Real>>({frequencies.size()});
  for (size_t i{0}; i < frequencies.size(); ++i) {
    std::complex<Real> sum{0.0, 0.0};

    for (size_t j{0}; j < time_domain_data.size(); ++j) {
      auto t{(j + 1) * dt + time_shift};
      sum += time_domain_data(j) *
             std::exp(-2.0 * constant::II * constant::PI * frequencies(i) * t);
    }
    res(i) = sum * dt;
  }
  return res;
}
}  // namespace xfdtd

#endif  // _XFDTD_CORE_DFT_H_
