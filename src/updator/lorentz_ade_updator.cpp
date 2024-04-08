#include <xfdtd/util/fdtd_basic.h>

#include "updator/dispersive_material_updator.h"

namespace xfdtd {

LorentzADEUpdator::LorentzADEUpdator(
    std::shared_ptr<const GridSpace> grid_space,
    std::shared_ptr<const CalculationParam> calculation_param,
    std::shared_ptr<EMF> emf, IndexTask task)
    : LinearDispersiveMaterialADEUpdator{std::move(grid_space),
                                         std::move(calculation_param),
                                         std::move(emf), task} {
  init();
}

void LorentzADEUpdator::init() {
  _lorentz_map = {};
  _lorentz_mediums = {};
  handleDispersiveMaterialADEUpdator(
      _lorentz_map, _lorentz_mediums,
      _calculation_param->materialParam()->materialArray());

  const auto& nx = task().xRange().size();
  const auto& ny = task().yRange().size();
  const auto& nz = task().zRange().size();
  const auto& nm = _lorentz_mediums.size();
  _coeff.resize(nm);
  _j_record.resize(nm);
  _ex_prev = xt::zeros<double>({nx, ny + 1, nz + 1});
  _ey_prev = xt::zeros<double>({nx + 1, ny, nz + 1});
  _ez_prev = xt::zeros<double>({nx + 1, ny + 1, nz});
  for (std::size_t i{0}; i < nm; ++i) {
    _coeff[i] = _lorentz_mediums[i]->coeffForADE();
    auto num_p = _lorentz_mediums[i]->numberOfPoles();
    _j_record[i]._jx = xt::zeros<double>({nx, ny + 1, nz + 1, num_p});
    _j_record[i]._jy = xt::zeros<double>({nx + 1, ny, nz + 1, num_p});
    _j_record[i]._jz = xt::zeros<double>({nx + 1, ny + 1, nz, num_p});
    _j_record[i]._jx_prev = xt::zeros<double>({nx, ny + 1, nz + 1, num_p});
    _j_record[i]._jy_prev = xt::zeros<double>({
        nx + 1,
        ny,
        nz + 1,
        num_p,
    });
    _j_record[i]._jz_prev = xt::zeros<double>({nx + 1, ny + 1, nz, num_p});
  }
}

void LorentzADEUpdator::updateE() {
  const auto task = this->task();
  const auto& cexe{_calculation_param->fdtdCoefficient()->cexe()};
  const auto& cexhy{_calculation_param->fdtdCoefficient()->cexhy()};
  const auto& cexhz{_calculation_param->fdtdCoefficient()->cexhz()};
  const auto& ceye{_calculation_param->fdtdCoefficient()->ceye()};
  const auto& ceyhz{_calculation_param->fdtdCoefficient()->ceyhz()};
  const auto& ceyhx{_calculation_param->fdtdCoefficient()->ceyhx()};
  const auto& ceze{_calculation_param->fdtdCoefficient()->ceze()};
  const auto& cezhx{_calculation_param->fdtdCoefficient()->cezhx()};
  const auto& cezhy{_calculation_param->fdtdCoefficient()->cezhy()};

  const auto& hx{_emf->hx()};
  const auto& hy{_emf->hy()};
  const auto& hz{_emf->hz()};
  auto& ex{_emf->ex()};
  auto& ey{_emf->ey()};
  auto& ez{_emf->ez()};

  auto& ex_prev = _ex_prev;
  auto& ey_prev = _ey_prev;
  auto& ez_prev = _ez_prev;

  const auto& dt = _calculation_param->timeParam()->dt();

  auto get_j_sum = [](std::size_t i, std::size_t j, std::size_t k,
                      const auto& coeff, const auto& j_arr,
                      const auto& j_prev_arr) {
    double j_sum{0};
    auto num_p = coeff._alpha.size();
    for (decltype(num_p) p = 0; p < num_p; ++p) {
      const auto& alpha = coeff._alpha(p);
      const auto& xi = coeff._xi(p);
      const auto& jj = j_arr(i, j, k, p);
      const auto& j_prev = j_prev_arr(i, j, k, p);
      j_sum += ((1 + alpha)) * jj + xi * j_prev;
    }
    return 0.5 * j_sum;
  };

  auto update_j = [](std::size_t i, std::size_t j, std::size_t k,
                     const auto& coeff, const auto& dt, const auto& e_next,
                     const auto& e_prev, auto&& jj, auto&& j_prev) {
    auto num_p = coeff._alpha.size();
    for (decltype(num_p) p = 0; p < num_p; ++p) {
      const auto& alpha = coeff._alpha(p);
      const auto& xi = coeff._xi(p);
      const auto& gamma = coeff._gamma(p);

      auto j_temp = jj(i, j, k, p);
      jj(i, j, k, p) = alpha * j_temp + xi * j_prev(i, j, k, p) +
                       gamma * (e_next - e_prev) / (2 * dt);
      j_prev(i, j, k, p) = j_temp;
    }
  };

  auto is = basic::GridStructure::exFDTDUpdateXStart(task.xRange().start());
  auto ie = basic::GridStructure::exFDTDUpdateXEnd(task.xRange().end());
  auto js = basic::GridStructure::exFDTDUpdateYStart(task.yRange().start());
  auto je = basic::GridStructure::exFDTDUpdateYEnd(task.yRange().end());
  auto ks = basic::GridStructure::exFDTDUpdateZStart(task.zRange().start());
  auto ke = basic::GridStructure::exFDTDUpdateZEnd(task.zRange().end());
  for (std::size_t i{is}; i < ie; ++i) {
    for (std::size_t j{js}; j < je; ++j) {
      for (std::size_t k{ks}; k < ke; ++k) {
        auto m_index =
            _grid_space->gridWithMaterial()(i, j, k)->materialIndex();
        if (m_index == -1 || _lorentz_map.find(m_index) == _lorentz_map.end()) {
          ex(i, j, k) = cexe(i, j, k) * ex(i, j, k) +
                        cexhy(i, j, k) * (hy(i, j, k) - hy(i, j, k - 1)) +
                        cexhz(i, j, k) * (hz(i, j, k) - hz(i, j - 1, k));
        } else {
          auto l_index = _lorentz_map.find(m_index)->second;
          const auto& coeff_ade = _coeff[l_index];
          auto& j_record = _j_record[l_index];

          double j_sum =
              get_j_sum(i, j, k, coeff_ade, j_record._jx, j_record._jx_prev);

          const auto& c1 = coeff_ade._c1;
          const auto& c2 = coeff_ade._c2;
          const auto& c3 = coeff_ade._c3;
          const auto e_prev_temp = ex(i, j, k);
          ex(i, j, k) = c1 * ex_prev(i, j, k) + cexe(i, j, k) * ex(i, j, k) +
                        cexhy(i, j, k) * (hy(i, j, k) - hy(i, j, k - 1)) +
                        cexhz(i, j, k) * (hz(i, j, k) - hz(i, j - 1, k)) -
                        c3 * j_sum;

          update_j(i, j, k, coeff_ade, dt, ex(i, j, k), ex_prev(i, j, k),
                   j_record._jx, j_record._jx_prev);
          ex_prev(i, j, k) = e_prev_temp;
        }
      }
    }
  }

  is = basic::GridStructure::eyFDTDUpdateXStart(task.xRange().start());
  ie = basic::GridStructure::eyFDTDUpdateXEnd(task.xRange().end());
  js = basic::GridStructure::eyFDTDUpdateYStart(task.yRange().start());
  je = basic::GridStructure::eyFDTDUpdateYEnd(task.yRange().end());
  ks = basic::GridStructure::eyFDTDUpdateZStart(task.zRange().start());
  ke = basic::GridStructure::eyFDTDUpdateZEnd(task.zRange().end());
  for (std::size_t i{is}; i < ie; ++i) {
    for (std::size_t j{js}; j < je; ++j) {
      for (std::size_t k{ks}; k < ke; ++k) {
        auto m_index =
            _grid_space->gridWithMaterial()(i, j, k)->materialIndex();

        if (m_index == -1 || _lorentz_map.find(m_index) == _lorentz_map.end()) {
          ey(i, j, k) = ceye(i, j, k) * ey(i, j, k) +
                        ceyhz(i, j, k) * (hz(i, j, k) - hz(i - 1, j, k)) +
                        ceyhx(i, j, k) * (hx(i, j, k) - hx(i, j, k - 1));
        } else {
          auto l_index = _lorentz_map.find(m_index)->second;
          const auto& coeff_ade = _coeff[l_index];
          auto& j_record = _j_record[l_index];

          double j_sum =
              get_j_sum(i, j, k, coeff_ade, j_record._jy, j_record._jy_prev);

          const auto& c1 = coeff_ade._c1;
          const auto& c2 = coeff_ade._c2;
          const auto& c3 = coeff_ade._c3;
          const auto e_prev_temp = ey(i, j, k);

          ey(i, j, k) = c1 * ey_prev(i, j, k) + ceye(i, j, k) * ey(i, j, k) +
                        ceyhz(i, j, k) * (hz(i, j, k) - hz(i - 1, j, k)) +
                        ceyhx(i, j, k) * (hx(i, j, k) - hx(i, j, k - 1)) -
                        c3 * j_sum;

          update_j(i, j, k, coeff_ade, dt, ey(i, j, k), ey_prev(i, j, k),
                   j_record._jy, j_record._jy_prev);
          ey_prev(i, j, k) = e_prev_temp;
        }
      }
    }
  }

  is = basic::GridStructure::ezFDTDUpdateXStart(task.xRange().start());
  ie = basic::GridStructure::ezFDTDUpdateXEnd(task.xRange().end());
  js = basic::GridStructure::ezFDTDUpdateYStart(task.yRange().start());
  je = basic::GridStructure::ezFDTDUpdateYEnd(task.yRange().end());
  ks = basic::GridStructure::ezFDTDUpdateZStart(task.zRange().start());
  ke = basic::GridStructure::ezFDTDUpdateZEnd(task.zRange().end());
  for (std::size_t i{is}; i < ie; ++i) {
    for (std::size_t j{js}; j < je; ++j) {
      for (std::size_t k{ks}; k < ke; ++k) {
        auto m_index =
            _grid_space->gridWithMaterial()(i, j, k)->materialIndex();

        if (m_index == -1 || _lorentz_map.find(m_index) == _lorentz_map.end()) {
          ez(i, j, k) = ceze(i, j, k) * ez(i, j, k) +
                        cezhx(i, j, k) * (hx(i, j, k) - hx(i, j - 1, k)) +
                        cezhy(i, j, k) * (hy(i, j, k) - hy(i - 1, j, k));
        } else {
          auto l_index = _lorentz_map.find(m_index)->second;
          const auto& coeff_ade = _coeff[l_index];
          auto& j_record = _j_record[l_index];

          double j_sum =
              get_j_sum(i, j, k, coeff_ade, j_record._jz, j_record._jz_prev);

          const auto& c1 = coeff_ade._c1;
          const auto& c2 = coeff_ade._c2;
          const auto& c3 = coeff_ade._c3;
          const auto e_prev_temp = ez(i, j, k);

          ez(i, j, k) = c1 * ez_prev(i, j, k) + ceze(i, j, k) * ez(i, j, k) +
                        cezhx(i, j, k) * (hx(i, j, k) - hx(i, j - 1, k)) +
                        cezhy(i, j, k) * (hy(i, j, k) - hy(i - 1, j, k)) -
                        c3 * j_sum;

          update_j(i, j, k, coeff_ade, dt, ez(i, j, k), ez_prev(i, j, k),
                   j_record._jz, j_record._jz_prev);
          ez_prev(i, j, k) = e_prev_temp;
        }
      }
    }
  }

  updateEEdge();
}

auto LorentzADEUpdator::updateEEdge() -> void {
  const auto is = task().xRange().start();
  const auto ie = task().xRange().end();
  const auto js = task().yRange().start();
  const auto je = task().yRange().end();
  const auto ks = task().zRange().start();
  const auto ke = task().zRange().end();

  const auto& cexe{_calculation_param->fdtdCoefficient()->cexe()};
  const auto& cexhy{_calculation_param->fdtdCoefficient()->cexhy()};
  const auto& cexhz{_calculation_param->fdtdCoefficient()->cexhz()};
  const auto& ceye{_calculation_param->fdtdCoefficient()->ceye()};
  const auto& ceyhz{_calculation_param->fdtdCoefficient()->ceyhz()};
  const auto& ceyhx{_calculation_param->fdtdCoefficient()->ceyhx()};
  const auto& ceze{_calculation_param->fdtdCoefficient()->ceze()};
  const auto& cezhx{_calculation_param->fdtdCoefficient()->cezhx()};
  const auto& cezhy{_calculation_param->fdtdCoefficient()->cezhy()};

  const auto& hx{_emf->hx()};
  const auto& hy{_emf->hy()};
  const auto& hz{_emf->hz()};
  auto& ex{_emf->ex()};
  auto& ey{_emf->ey()};
  auto& ez{_emf->ez()};

  auto& ex_prev = _ex_prev;
  auto& ey_prev = _ey_prev;
  auto& ez_prev = _ez_prev;

  const auto& dt = _calculation_param->timeParam()->dt();

  bool contain_xn_edge = containXNEdge();
  bool contain_yn_edge = containYNEdge();
  bool contain_zn_edge = containZNEdge();

  auto get_j_sum = [](std::size_t i, std::size_t j, std::size_t k,
                      const auto& coeff, const auto& j_arr,
                      const auto& j_prev_arr) {
    double j_sum{0};
    auto num_p = coeff._alpha.size();
    for (decltype(num_p) p = 0; p < num_p; ++p) {
      const auto& alpha = coeff._alpha(p);
      const auto& xi = coeff._xi(p);
      const auto& jj = j_arr(i, j, k, p);
      const auto& j_prev = j_prev_arr(i, j, k, p);
      j_sum += ((1 + alpha)) * jj + xi * j_prev;
    }
    return 0.5 * j_sum;
  };

  auto update_j = [](std::size_t i, std::size_t j, std::size_t k,
                     const auto& coeff, const auto& dt, const auto& e_next,
                     const auto& e_prev, auto&& jj, auto&& j_prev) {
    auto num_p = coeff._alpha.size();
    for (decltype(num_p) p = 0; p < num_p; ++p) {
      const auto& alpha = coeff._alpha(p);
      const auto& xi = coeff._xi(p);
      const auto& gamma = coeff._gamma(p);

      auto j_temp = jj(i, j, k, p);
      jj(i, j, k, p) = alpha * j_temp + xi * j_prev(i, j, k, p) +
                       gamma * (e_next - e_prev) / (2 * dt);
      j_prev(i, j, k, p) = j_temp;
    }
  };

  if (!contain_yn_edge && !contain_zn_edge) {
    auto j = js;
    auto k = ks;
    for (std::size_t i{is}; i < ie; ++i) {
      auto m_index = _grid_space->gridWithMaterial()(i, j, k)->materialIndex();
      if (m_index == -1 || _lorentz_map.find(m_index) == _lorentz_map.end()) {
        ex(i, j, k) = cexe(i, j, k) * ex(i, j, k) +
                      cexhy(i, j, k) * (hy(i, j, k) - hy(i, j, k - 1)) +
                      cexhz(i, j, k) * (hz(i, j, k) - hz(i, j - 1, k));
      } else {
        auto l_index = _lorentz_map.find(m_index)->second;
        const auto& coeff_ade = _coeff[l_index];
        auto& j_record = _j_record[l_index];

        double j_sum =
            get_j_sum(i, j, k, coeff_ade, j_record._jx, j_record._jx_prev);

        const auto& c1 = coeff_ade._c1;
        const auto& c2 = coeff_ade._c2;
        const auto& c3 = coeff_ade._c3;
        const auto e_prev_temp = ex(i, j, k);
        ex(i, j, k) = c1 * ex_prev(i, j, k) + cexe(i, j, k) * ex(i, j, k) +
                      cexhy(i, j, k) * (hy(i, j, k) - hy(i, j, k - 1)) +
                      cexhz(i, j, k) * (hz(i, j, k) - hz(i, j - 1, k)) -
                      c3 * j_sum;

        update_j(i, j, k, coeff_ade, dt, ex(i, j, k), ex_prev(i, j, k),
                 j_record._jx, j_record._jx_prev);
        ex_prev(i, j, k) = e_prev_temp;
      }
    }
  }

  if (!contain_xn_edge && !contain_zn_edge) {
    auto i = is;
    auto k = ks;
    for (std::size_t j{js}; j < je; ++j) {
      auto m_index = _grid_space->gridWithMaterial()(i, j, k)->materialIndex();

      if (m_index == -1 || _lorentz_map.find(m_index) == _lorentz_map.end()) {
        ey(i, j, k) = ceye(i, j, k) * ey(i, j, k) +
                      ceyhz(i, j, k) * (hz(i, j, k) - hz(i - 1, j, k)) +
                      ceyhx(i, j, k) * (hx(i, j, k) - hx(i, j, k - 1));
      } else {
        auto l_index = _lorentz_map.find(m_index)->second;
        const auto& coeff_ade = _coeff[l_index];
        auto& j_record = _j_record[l_index];

        double j_sum =
            get_j_sum(i, j, k, coeff_ade, j_record._jy, j_record._jy_prev);

        const auto& c1 = coeff_ade._c1;
        const auto& c2 = coeff_ade._c2;
        const auto& c3 = coeff_ade._c3;
        const auto e_prev_temp = ey(i, j, k);

        ey(i, j, k) = c1 * ey_prev(i, j, k) + ceye(i, j, k) * ey(i, j, k) +
                      ceyhz(i, j, k) * (hz(i, j, k) - hz(i - 1, j, k)) +
                      ceyhx(i, j, k) * (hx(i, j, k) - hx(i, j, k - 1)) -
                      c3 * j_sum;

        update_j(i, j, k, coeff_ade, dt, ey(i, j, k), ey_prev(i, j, k),
                 j_record._jy, j_record._jy_prev);
        ey_prev(i, j, k) = e_prev_temp;
      }
    }
  }

  if (!contain_xn_edge && !contain_yn_edge) {
    auto i = is;
    auto j = js;
    for (std::size_t k{ks}; k < ke; ++k) {
      auto m_index = _grid_space->gridWithMaterial()(i, j, k)->materialIndex();

      if (m_index == -1 || _lorentz_map.find(m_index) == _lorentz_map.end()) {
        ez(i, j, k) = ceze(i, j, k) * ez(i, j, k) +
                      cezhx(i, j, k) * (hx(i, j, k) - hx(i, j - 1, k)) +
                      cezhy(i, j, k) * (hy(i, j, k) - hy(i - 1, j, k));
      } else {
        auto l_index = _lorentz_map.find(m_index)->second;
        const auto& coeff_ade = _coeff[l_index];
        auto& j_record = _j_record[l_index];

        double j_sum =
            get_j_sum(i, j, k, coeff_ade, j_record._jz, j_record._jz_prev);

        const auto& c1 = coeff_ade._c1;
        const auto& c2 = coeff_ade._c2;
        const auto& c3 = coeff_ade._c3;
        const auto e_prev_temp = ez(i, j, k);

        ez(i, j, k) = c1 * ez_prev(i, j, k) + ceze(i, j, k) * ez(i, j, k) +
                      cezhx(i, j, k) * (hx(i, j, k) - hx(i, j - 1, k)) +
                      cezhy(i, j, k) * (hy(i, j, k) - hy(i - 1, j, k)) -
                      c3 * j_sum;

        update_j(i, j, k, coeff_ade, dt, ez(i, j, k), ez_prev(i, j, k),
                 j_record._jz, j_record._jz_prev);
        ez_prev(i, j, k) = e_prev_temp;
      }
    }
  }

  if (!contain_xn_edge) {
    auto i = is;
    for (std::size_t j{js}; j < je; ++j) {
      for (std::size_t k{ks + 1}; k < ke; ++k) {
        auto m_index =
            _grid_space->gridWithMaterial()(i, j, k)->materialIndex();

        if (m_index == -1 || _lorentz_map.find(m_index) == _lorentz_map.end()) {
          ey(i, j, k) = ceye(i, j, k) * ey(i, j, k) +
                        ceyhz(i, j, k) * (hz(i, j, k) - hz(i - 1, j, k)) +
                        ceyhx(i, j, k) * (hx(i, j, k) - hx(i, j, k - 1));
        } else {
          auto l_index = _lorentz_map.find(m_index)->second;
          const auto& coeff_ade = _coeff[l_index];
          auto& j_record = _j_record[l_index];

          double j_sum =
              get_j_sum(i, j, k, coeff_ade, j_record._jy, j_record._jy_prev);

          const auto& c1 = coeff_ade._c1;
          const auto& c2 = coeff_ade._c2;
          const auto& c3 = coeff_ade._c3;
          const auto e_prev_temp = ey(i, j, k);

          ey(i, j, k) = c1 * ey_prev(i, j, k) + ceye(i, j, k) * ey(i, j, k) +
                        ceyhz(i, j, k) * (hz(i, j, k) - hz(i - 1, j, k)) +
                        ceyhx(i, j, k) * (hx(i, j, k) - hx(i, j, k - 1)) -
                        c3 * j_sum;

          update_j(i, j, k, coeff_ade, dt, ey(i, j, k), ey_prev(i, j, k),
                   j_record._jy, j_record._jy_prev);
          ey_prev(i, j, k) = e_prev_temp;
        }
      }
    }
    for (std::size_t j{js + 1}; j < je; ++j) {
      for (std::size_t k{ks}; k < ke; ++k) {
        auto m_index =
            _grid_space->gridWithMaterial()(i, j, k)->materialIndex();

        if (m_index == -1 || _lorentz_map.find(m_index) == _lorentz_map.end()) {
          ez(i, j, k) = ceze(i, j, k) * ez(i, j, k) +
                        cezhx(i, j, k) * (hx(i, j, k) - hx(i, j - 1, k)) +
                        cezhy(i, j, k) * (hy(i, j, k) - hy(i - 1, j, k));
        } else {
          auto l_index = _lorentz_map.find(m_index)->second;
          const auto& coeff_ade = _coeff[l_index];
          auto& j_record = _j_record[l_index];

          double j_sum =
              get_j_sum(i, j, k, coeff_ade, j_record._jz, j_record._jz_prev);

          const auto& c1 = coeff_ade._c1;
          const auto& c2 = coeff_ade._c2;
          const auto& c3 = coeff_ade._c3;
          const auto e_prev_temp = ez(i, j, k);

          ez(i, j, k) = c1 * ez_prev(i, j, k) + ceze(i, j, k) * ez(i, j, k) +
                        cezhx(i, j, k) * (hx(i, j, k) - hx(i, j - 1, k)) +
                        cezhy(i, j, k) * (hy(i, j, k) - hy(i - 1, j, k)) -
                        c3 * j_sum;

          update_j(i, j, k, coeff_ade, dt, ez(i, j, k), ez_prev(i, j, k),
                   j_record._jz, j_record._jz_prev);
          ez_prev(i, j, k) = e_prev_temp;
        }
      }
    }
  }

  if (!contain_yn_edge) {
    auto j = js;
    for (std::size_t i{is + 1}; i < ie; ++i) {
      for (std::size_t k{ks}; k < ke; ++k) {
        auto m_index =
            _grid_space->gridWithMaterial()(i, j, k)->materialIndex();

        if (m_index == -1 || _lorentz_map.find(m_index) == _lorentz_map.end()) {
          ez(i, j, k) = ceze(i, j, k) * ez(i, j, k) +
                        cezhx(i, j, k) * (hx(i, j, k) - hx(i, j - 1, k)) +
                        cezhy(i, j, k) * (hy(i, j, k) - hy(i - 1, j, k));
        } else {
          auto l_index = _lorentz_map.find(m_index)->second;
          const auto& coeff_ade = _coeff[l_index];
          auto& j_record = _j_record[l_index];

          double j_sum =
              get_j_sum(i, j, k, coeff_ade, j_record._jz, j_record._jz_prev);

          const auto& c1 = coeff_ade._c1;
          const auto& c2 = coeff_ade._c2;
          const auto& c3 = coeff_ade._c3;
          const auto e_prev_temp = ez(i, j, k);

          ez(i, j, k) = c1 * ez_prev(i, j, k) + ceze(i, j, k) * ez(i, j, k) +
                        cezhx(i, j, k) * (hx(i, j, k) - hx(i, j - 1, k)) +
                        cezhy(i, j, k) * (hy(i, j, k) - hy(i - 1, j, k)) -
                        c3 * j_sum;

          update_j(i, j, k, coeff_ade, dt, ez(i, j, k), ez_prev(i, j, k),
                   j_record._jz, j_record._jz_prev);
          ez_prev(i, j, k) = e_prev_temp;
        }
      }
    }
    for (std::size_t i{is}; i < ie; ++i) {
      for (std::size_t k{ks + 1}; k < ke; ++k) {
        auto m_index =
            _grid_space->gridWithMaterial()(i, j, k)->materialIndex();
        if (m_index == -1 || _lorentz_map.find(m_index) == _lorentz_map.end()) {
          ex(i, j, k) = cexe(i, j, k) * ex(i, j, k) +
                        cexhy(i, j, k) * (hy(i, j, k) - hy(i, j, k - 1)) +
                        cexhz(i, j, k) * (hz(i, j, k) - hz(i, j - 1, k));
        } else {
          auto l_index = _lorentz_map.find(m_index)->second;
          const auto& coeff_ade = _coeff[l_index];
          auto& j_record = _j_record[l_index];

          double j_sum =
              get_j_sum(i, j, k, coeff_ade, j_record._jx, j_record._jx_prev);

          const auto& c1 = coeff_ade._c1;
          const auto& c2 = coeff_ade._c2;
          const auto& c3 = coeff_ade._c3;
          const auto e_prev_temp = ex(i, j, k);
          ex(i, j, k) = c1 * ex_prev(i, j, k) + cexe(i, j, k) * ex(i, j, k) +
                        cexhy(i, j, k) * (hy(i, j, k) - hy(i, j, k - 1)) +
                        cexhz(i, j, k) * (hz(i, j, k) - hz(i, j - 1, k)) -
                        c3 * j_sum;

          update_j(i, j, k, coeff_ade, dt, ex(i, j, k), ex_prev(i, j, k),
                   j_record._jx, j_record._jx_prev);
          ex_prev(i, j, k) = e_prev_temp;
        }
      }
    }
  }

  if (!contain_zn_edge) {
    auto k = ks;
    for (std::size_t i{is}; i < ie; ++i) {
      for (std::size_t j{js + 1}; j < je; ++j) {
        auto m_index =
            _grid_space->gridWithMaterial()(i, j, k)->materialIndex();
        if (m_index == -1 || _lorentz_map.find(m_index) == _lorentz_map.end()) {
          ex(i, j, k) = cexe(i, j, k) * ex(i, j, k) +
                        cexhy(i, j, k) * (hy(i, j, k) - hy(i, j, k - 1)) +
                        cexhz(i, j, k) * (hz(i, j, k) - hz(i, j - 1, k));
        } else {
          auto l_index = _lorentz_map.find(m_index)->second;
          const auto& coeff_ade = _coeff[l_index];
          auto& j_record = _j_record[l_index];

          double j_sum =
              get_j_sum(i, j, k, coeff_ade, j_record._jx, j_record._jx_prev);

          const auto& c1 = coeff_ade._c1;
          const auto& c2 = coeff_ade._c2;
          const auto& c3 = coeff_ade._c3;
          const auto e_prev_temp = ex(i, j, k);
          ex(i, j, k) = c1 * ex_prev(i, j, k) + cexe(i, j, k) * ex(i, j, k) +
                        cexhy(i, j, k) * (hy(i, j, k) - hy(i, j, k - 1)) +
                        cexhz(i, j, k) * (hz(i, j, k) - hz(i, j - 1, k)) -
                        c3 * j_sum;

          update_j(i, j, k, coeff_ade, dt, ex(i, j, k), ex_prev(i, j, k),
                   j_record._jx, j_record._jx_prev);
          ex_prev(i, j, k) = e_prev_temp;
        }
      }
    }
    for (std::size_t i{is + 1}; i < ie; ++i) {
      for (std::size_t j{js}; j < je; ++j) {
        auto m_index =
            _grid_space->gridWithMaterial()(i, j, k)->materialIndex();

        if (m_index == -1 || _lorentz_map.find(m_index) == _lorentz_map.end()) {
          ey(i, j, k) = ceye(i, j, k) * ey(i, j, k) +
                        ceyhz(i, j, k) * (hz(i, j, k) - hz(i - 1, j, k)) +
                        ceyhx(i, j, k) * (hx(i, j, k) - hx(i, j, k - 1));
        } else {
          auto l_index = _lorentz_map.find(m_index)->second;
          const auto& coeff_ade = _coeff[l_index];
          auto& j_record = _j_record[l_index];

          double j_sum =
              get_j_sum(i, j, k, coeff_ade, j_record._jy, j_record._jy_prev);

          const auto& c1 = coeff_ade._c1;
          const auto& c2 = coeff_ade._c2;
          const auto& c3 = coeff_ade._c3;
          const auto e_prev_temp = ey(i, j, k);

          ey(i, j, k) = c1 * ey_prev(i, j, k) + ceye(i, j, k) * ey(i, j, k) +
                        ceyhz(i, j, k) * (hz(i, j, k) - hz(i - 1, j, k)) +
                        ceyhx(i, j, k) * (hx(i, j, k) - hx(i, j, k - 1)) -
                        c3 * j_sum;

          update_j(i, j, k, coeff_ade, dt, ey(i, j, k), ey_prev(i, j, k),
                   j_record._jy, j_record._jy_prev);
          ey_prev(i, j, k) = e_prev_temp;
        }
      }
    }
  }
}

}  // namespace xfdtd
