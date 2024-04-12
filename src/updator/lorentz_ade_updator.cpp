#include <xfdtd/util/fdtd_basic.h>

#include "updator/dispersive_material_updator.h"
#include "updator/update_scheme.h"

namespace xfdtd {

LorentzADECorrector::LorentzADECorrector(
    const LorentzMedium& lorentz_medium,
    std::shared_ptr<const CalculationParam> calculation_param,
    std::shared_ptr<EMF> emf, IndexTask task)
    : ADECorrector{lorentz_medium.numberOfPoles(), task,
                   lorentz_medium.coeffForADE()._c3,
                   std::move(calculation_param), std::move(emf)},
      _jx_prev{xt::zeros<Real>({task.xRange().size(), task.yRange().size() + 1,
                                task.zRange().size() + 1, numPole()})},
      _jy_prev{xt::zeros<Real>({task.xRange().size() + 1, task.yRange().size(),
                                task.zRange().size() + 1, numPole()})},
      _jz_prev{
          xt::zeros<Real>({task.xRange().size() + 1, task.yRange().size() + 1,
                           task.zRange().size(), numPole()})},
      _ex_prev{xt::zeros<Real>({task.xRange().size(), task.yRange().size() + 1,
                                task.zRange().size() + 1})},
      _ey_prev{xt::zeros<Real>({task.xRange().size() + 1, task.yRange().size(),
                                task.zRange().size() + 1})},
      _ez_prev{
          xt::zeros<Real>({task.xRange().size() + 1, task.yRange().size() + 1,
                           task.zRange().size()})},
      _alpha{lorentz_medium.coeffForADE()._alpha},
      _xi{lorentz_medium.coeffForADE()._xi},
      _gamma{lorentz_medium.coeffForADE()._gamma},
      _c1{lorentz_medium.coeffForADE()._c1} {}

auto LorentzADECorrector::updateEx(Index node_i, Index node_j, Index node_k)
    -> void {
  auto& ex{_emf->ex()};
  const auto& cexe{_calculation_param->fdtdCoefficient()->cexe()};
  const auto& cexhy{_calculation_param->fdtdCoefficient()->cexhy()};
  const auto& cexhz{_calculation_param->fdtdCoefficient()->cexhz()};
  const auto& hy{_emf->hy()};
  const auto& hz{_emf->hz()};
  const auto coeff_j = coeffJ();
  const auto dt{_calculation_param->timeParam()->dt()};

  auto local_i = node_i - task().xRange().start();
  auto local_j = node_j - task().yRange().start();
  auto local_k = node_k - task().zRange().start();

  auto j_sum = calculateJSum(node_i, node_j, node_k, _jx, _jx_prev);
  const auto e_cur = ex(node_i, node_j, node_k);
  auto& e_next = ex(node_i, node_j, node_k);
  auto& e_prev = _ex_prev(local_i, local_j, local_k);

  e_next = _c1 * e_prev +
           eNext(cexe(node_i, node_j, node_k), ex(node_i, node_j, node_k),
                 cexhy(node_i, node_j, node_k), hy(node_i, node_j, node_k),
                 hy(node_i, node_j, node_k - 1), cexhz(node_i, node_j, node_k),
                 hz(node_i, node_j, node_k), hz(node_i, node_j - 1, node_k)) -
           coeff_j * j_sum;

  updateJ(node_i, node_j, node_k, e_next, e_prev, dt, _jx, _jx_prev);
  e_prev = e_cur;
}

auto LorentzADECorrector::updateEy(Index node_i, Index node_j, Index node_k)
    -> void {
  auto& ey{_emf->ey()};
  const auto& ceye{_calculation_param->fdtdCoefficient()->ceye()};
  const auto& ceyhz{_calculation_param->fdtdCoefficient()->ceyhz()};
  const auto& ceyhx{_calculation_param->fdtdCoefficient()->ceyhx()};
  const auto& hz{_emf->hz()};
  const auto& hx{_emf->hx()};
  const auto dt{_calculation_param->timeParam()->dt()};

  auto local_i = node_i - task().xRange().start();
  auto local_j = node_j - task().yRange().start();
  auto local_k = node_k - task().zRange().start();

  auto j_sum = calculateJSum(node_i, node_j, node_k, _jy, _jy_prev);
  const auto e_cur = ey(node_i, node_j, node_k);
  auto& e_next = ey(node_i, node_j, node_k);
  auto& e_prev = _ey_prev(local_i, local_j, local_k);

  e_next = _c1 * e_prev +
           eNext(ceye(node_i, node_j, node_k), ey(node_i, node_j, node_k),
                 ceyhz(node_i, node_j, node_k), hz(node_i, node_j, node_k),
                 hz(node_i - 1, node_j, node_k), ceyhx(node_i, node_j, node_k),
                 hx(node_i, node_j, node_k), hx(node_i, node_j, node_k - 1)) -
           coeffJ() * j_sum;

  updateJ(node_i, node_j, node_k, e_next, e_prev, dt, _jy, _jy_prev);
  e_prev = e_cur;
}

auto LorentzADECorrector::updateEz(Index node_i, Index node_j, Index node_k)
    -> void {
  auto& ez{_emf->ez()};
  const auto& ceze{_calculation_param->fdtdCoefficient()->ceze()};
  const auto& cezhx{_calculation_param->fdtdCoefficient()->cezhx()};
  const auto& cezhy{_calculation_param->fdtdCoefficient()->cezhy()};
  const auto& hx{_emf->hx()};
  const auto& hy{_emf->hy()};
  const auto dt{_calculation_param->timeParam()->dt()};

  auto local_i = node_i - task().xRange().start();
  auto local_j = node_j - task().yRange().start();
  auto local_k = node_k - task().zRange().start();

  auto j_sum = calculateJSum(node_i, node_j, node_k, _jz, _jz_prev);
  const auto e_cur = ez(node_i, node_j, node_k);
  auto& e_next = ez(node_i, node_j, node_k);
  auto& e_prev = _ez_prev(local_i, local_j, local_k);

  e_next = _c1 * e_prev +
           eNext(ceze(node_i, node_j, node_k), ez(node_i, node_j, node_k),
                 cezhx(node_i, node_j, node_k), hx(node_i, node_j, node_k),
                 hx(node_i, node_j - 1, node_k), cezhy(node_i, node_j, node_k),
                 hy(node_i, node_j, node_k), hy(node_i - 1, node_j, node_k)) -
           coeffJ() * j_sum;

  updateJ(node_i, node_j, node_k, e_next, e_prev, dt, _jz, _jz_prev);
  e_prev = e_cur;
}

auto LorentzADECorrector::calculateJSum(Index node_i, Index node_j,
                                        Index node_k,
                                        const Array4D<Real>& j_arr,
                                        const Array4D<Real>& j_prev_arr) const
    -> Real {
  Real j_sum{0};
  auto num_p = numPole();

  auto local_i = node_i - task().xRange().start();
  auto local_j = node_j - task().yRange().start();
  auto local_k = node_k - task().zRange().start();

  for (Index p{0}; p < num_p; ++p) {
    const auto alpha = _alpha(p);
    const auto xi = _xi(p);
    j_sum += ((1 + alpha)) * j_arr(local_i, local_j, local_k, p) +
             xi * j_prev_arr(local_i, local_j, local_k, p);
  }

  return 0.5 * j_sum;
}

auto LorentzADECorrector::updateJ(Index node_i, Index node_j, Index node_k,
                                  Real e_next, Real e_prev, Real dt,
                                  Array4D<Real>& j_arr,
                                  Array4D<Real>& j_prev_arr) -> void {
  auto num_p = numPole();

  auto local_i = node_i - task().xRange().start();
  auto local_j = node_j - task().yRange().start();
  auto local_k = node_k - task().zRange().start();

  for (Index p{0}; p < num_p; ++p) {
    const auto alpha = _alpha(p);
    const auto xi = _xi(p);
    const auto gamma = _gamma(p);

    const auto j_temp = j_arr(local_i, local_j, local_k, p);

    j_arr(local_i, local_j, local_k, p) =
        alpha * j_temp + xi * j_prev_arr(local_i, local_j, local_k, p) +
        gamma * (e_next - e_prev) / (2 * dt);
    j_prev_arr(local_i, local_j, local_k, p) = j_temp;
  }
}

}  // namespace xfdtd
