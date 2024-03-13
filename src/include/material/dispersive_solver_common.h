#ifndef _XFDTD_CORE_DISPERSIVE_SOLVER_COMMON_H_
#define _XFDTD_CORE_DISPERSIVE_SOLVER_COMMON_H_

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

struct DrudeCoeff {
  xt::xarray<double> _k;
  xt::xarray<double> _beta;

  double _a;
  double _b;
};

struct DebyCoeff {
  xt::xarray<double> _k;
  xt::xarray<double> _beta;

  double _a;
  double _b;
};

}  // namespace ade

namespace z_transfer {}

}  // namespace xfdtd

#endif  // _XFDTD_CORE_DISPERSIVE_SOLVER_COMMON_H_
