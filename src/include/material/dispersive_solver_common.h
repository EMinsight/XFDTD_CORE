#ifndef _XFDTD_LIB_DISPERSIVE_SOLVER_COMMON_H_
#define _XFDTD_LIB_DISPERSIVE_SOLVER_COMMON_H_

#include <xtensor/xarray.hpp>

namespace xfdtd {

namespace ade {
struct LorentzCoeff {
  xt::xarray<double> _alpha;
  xt::xarray<double> _xi;
  xt::xarray<double> _gamma;

  double _c1;
  double _c2;
  double _c3;
};

struct DrudeCoeff {};

}  // namespace ade

namespace z {}

}  // namespace xfdtd

#endif  // _XFDTD_LIB_DISPERSIVE_SOLVER_COMMON_H_
