#include <xfdtd/common/type_define.h>
#include <xfdtd/util/fdtd_basic.h>

#include "updator/dispersive_material_updator.h"
#include "updator/update_scheme.h"

namespace xfdtd {

DebyeADECorrector::DebyeADECorrector(
    const DebyeMedium& debye_medium,
    std::shared_ptr<const CalculationParam> calculation_param,
    std::shared_ptr<EMF> emf, IndexTask task)
    : LinearDispersiveMaterialADEUpdator::
          ADECorrector{debye_medium.numberOfPoles(), task,
                       debye_medium.coeffForADE()._b,
                       std::move(calculation_param), std::move(emf)},
      _k{debye_medium.coeffForADE()._k},
      _beta{debye_medium.coeffForADE()._beta} {}

auto DebyeADECorrector::updateEx(Index node_i, Index node_j, Index node_k)
    -> void {
  auto& ex{_emf->ex()};
  const auto& cexe{_calculation_param->fdtdCoefficient()->cexe()};
  const auto& cexhy{_calculation_param->fdtdCoefficient()->cexhy()};
  const auto& cexhz{_calculation_param->fdtdCoefficient()->cexhz()};
  const auto& hy{_emf->hy()};
  const auto& hz{_emf->hz()};
  const auto coeff_j = coeffJ();

  auto j_sum = calculateJSum(node_i, node_j, node_k, _jx);
  const auto e_cur = ex(node_i, node_j, node_k);
  auto& e_next = ex(node_i, node_j, node_k);

  ex(node_i, node_j, node_k) =
      eNext(cexe(node_i, node_j, node_k), ex(node_i, node_j, node_k),
            cexhy(node_i, node_j, node_k), hy(node_i, node_j, node_k),
            hy(node_i, node_j, node_k - 1), cexhz(node_i, node_j, node_k),
            hz(node_i, node_j, node_k), hz(node_i, node_j - 1, node_k)) -
      coeff_j * j_sum;

  const auto dt = _calculation_param->timeParam()->dt();
  updateJ(node_i, node_j, node_k, e_next, e_cur, dt, _jx);
}

auto DebyeADECorrector::updateEy(Index node_i, Index node_j, Index node_k)
    -> void {
  auto& ey{_emf->ey()};
  const auto& ceye{_calculation_param->fdtdCoefficient()->ceye()};
  const auto& ceyhz{_calculation_param->fdtdCoefficient()->ceyhz()};
  const auto& ceyhx{_calculation_param->fdtdCoefficient()->ceyhx()};
  const auto& hz{_emf->hz()};
  const auto& hx{_emf->hx()};

  auto j_sum = calculateJSum(node_i, node_j, node_k, _jy);
  const auto e_cur = ey(node_i, node_j, node_k);
  auto&& e_next = ey(node_i, node_j, node_k);

  e_next = eNext(ceye(node_i, node_j, node_k), ey(node_i, node_j, node_k),
                 ceyhz(node_i, node_j, node_k), hz(node_i, node_j, node_k),
                 hz(node_i - 1, node_j, node_k), ceyhx(node_i, node_j, node_k),
                 hx(node_i, node_j, node_k), hx(node_i, node_j, node_k - 1)) -
           coeffJ() * j_sum;

  const auto dt = _calculation_param->timeParam()->dt();
  updateJ(node_i, node_j, node_k, e_next, e_cur, dt, _jy);
}

auto DebyeADECorrector::updateEz(Index node_i, Index node_j, Index node_k)
    -> void {
  auto& ez{_emf->ez()};
  const auto& ceze{_calculation_param->fdtdCoefficient()->ceze()};
  const auto& cezhx{_calculation_param->fdtdCoefficient()->cezhx()};
  const auto& cezhy{_calculation_param->fdtdCoefficient()->cezhy()};
  const auto& hx{_emf->hx()};
  const auto& hy{_emf->hy()};

  auto j_sum = calculateJSum(node_i, node_j, node_k, _jz);
  const auto e_cur = ez(node_i, node_j, node_k);
  auto&& e_next = ez(node_i, node_j, node_k);

  e_next = eNext(ceze(node_i, node_j, node_k), ez(node_i, node_j, node_k),
                 cezhx(node_i, node_j, node_k), hx(node_i, node_j, node_k),
                 hx(node_i, node_j - 1, node_k), cezhy(node_i, node_j, node_k),
                 hy(node_i, node_j, node_k), hy(node_i - 1, node_j, node_k)) -
           coeffJ() * j_sum;

  const auto dt = _calculation_param->timeParam()->dt();
  updateJ(node_i, node_j, node_k, e_next, e_cur, dt, _jz);
}

auto DebyeADECorrector::calculateJSum(Index node_i, Index node_j, Index node_k,
                                      const Array4D<Real>& j_arr) const
    -> Real {
  Real j_sum{0};
  auto num_p = numPole();

  auto local_i = node_i - task().xRange().start();
  auto local_j = node_j - task().yRange().start();
  auto local_k = node_k - task().zRange().start();

  for (Index p{0}; p < num_p; ++p) {
    auto k_p = _k(p);
    j_sum += (1 + k_p) * j_arr(local_i, local_j, local_k, p);
  }

  return 0.5 * j_sum;
}

auto DebyeADECorrector::updateJ(Index node_i, Index node_j, Index node_k,
                                Real e_next, Real e_cur, Real dt,
                                Array4D<Real>& j_arr) -> void {
  auto num_p = numPole();
  auto local_i = node_i - task().xRange().start();
  auto local_j = node_j - task().yRange().start();
  auto local_k = node_k - task().zRange().start();

  for (Index p{0}; p < num_p; ++p) {
    const auto& k_p = _k(p);
    const auto& beta = _beta(p);

    j_arr(local_i, local_j, local_k, p) =
        k_p * j_arr(local_i, local_j, local_k, p) +
        beta * (e_next - e_cur) / dt;
  }
}

}  // namespace xfdtd
