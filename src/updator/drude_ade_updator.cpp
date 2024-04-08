#include <xfdtd/util/fdtd_basic.h>

#include "updator/dispersive_material_updator.h"

namespace xfdtd {

DrudeADEUpdator::DrudeADEUpdator(
    std::shared_ptr<const GridSpace> grid_space,
    std::shared_ptr<const CalculationParam> calculation_param,
    std::shared_ptr<EMF> emf, IndexTask task)
    : LinearDispersiveMaterialADEUpdator{std::move(grid_space),
                                         std::move(calculation_param),
                                         std::move(emf), task} {
  init();
}

void DrudeADEUpdator::init() {
  const auto& nx = task().xRange().size();
  const auto& ny = task().yRange().size();
  const auto& nz = task().zRange().size();

  _drude_map = {};
  _drude_mediums = {};

  handleDispersiveMaterialADEUpdator(
      _drude_map, _drude_mediums,
      _calculation_param->materialParam()->materialArray());

  const auto& nm = _drude_mediums.size();
  _coeff.resize(nm);
  _j_record.resize(nm);
  for (std::size_t i{0}; i < nm; ++i) {
    _coeff[i] = _drude_mediums[i]->coeffForADE();
    auto num_p = _drude_mediums[i]->numberOfPoles();
    _j_record[i]._jx = xt::zeros<double>({nx, ny + 1, nz + 1, num_p});
    _j_record[i]._jy = xt::zeros<double>({nx + 1, ny, nz + 1, num_p});
    _j_record[i]._jz = xt::zeros<double>({nx + 1, ny + 1, nz, num_p});
  }
}

void DrudeADEUpdator::updateE() {
  const auto task = this->task();

  auto& cexe{_calculation_param->fdtdCoefficient()->cexe()};
  auto& cexhy{_calculation_param->fdtdCoefficient()->cexhy()};
  auto& cexhz{_calculation_param->fdtdCoefficient()->cexhz()};
  auto& ceye{_calculation_param->fdtdCoefficient()->ceye()};
  auto& ceyhz{_calculation_param->fdtdCoefficient()->ceyhz()};
  auto& ceyhx{_calculation_param->fdtdCoefficient()->ceyhx()};
  auto& ceze{_calculation_param->fdtdCoefficient()->ceze()};
  auto& cezhx{_calculation_param->fdtdCoefficient()->cezhx()};
  auto& cezhy{_calculation_param->fdtdCoefficient()->cezhy()};

  auto& hx{_emf->hx()};
  auto& hy{_emf->hy()};
  auto& hz{_emf->hz()};
  auto& ex{_emf->ex()};
  auto& ey{_emf->ey()};
  auto& ez{_emf->ez()};

  const auto& dt = _calculation_param->timeParam()->dt();

  auto get_j_sum = [](std::size_t i, std::size_t j, std::size_t k,
                      const auto& coeff, const auto& j_arr) {
    double j_sum{0};
    auto num_p = coeff._k.size();
    for (decltype(num_p) p = 0; p < num_p; ++p) {
      const auto& k = coeff._k(p);
      const auto& jj = j_arr(i, j, k, p);
      j_sum += ((1 + k)) * jj;
    }
    return 0.5 * j_sum;
  };

  auto update_j = [](std::size_t i, std::size_t j, std::size_t k,
                     const auto& coeff, const auto& dt, const auto& e_next,
                     const auto& e_cur, auto&& jj) {
    auto num_p = coeff._k.size();
    for (decltype(num_p) p = 0; p < num_p; ++p) {
      const auto& k = coeff._k(p);
      const auto& beta = coeff._beta(p);

      jj(i, j, k, p) = k * jj(i, j, k, p) + beta * (e_next + e_cur);
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
        if (m_index == -1 || _drude_map.find(m_index) == _drude_map.end()) {
          ex(i, j, k) = cexe(i, j, k) * ex(i, j, k) +
                        cexhy(i, j, k) * (hy(i, j, k) - hy(i, j, k - 1)) +
                        cexhz(i, j, k) * (hz(i, j, k) - hz(i, j - 1, k));
        } else {
          auto l_index = _drude_map.find(m_index)->second;
          const auto& coeff_ade = _coeff[l_index];
          auto& j_record = _j_record[l_index];

          double j_sum = get_j_sum(i, j, k, coeff_ade, j_record._jx);

          const auto& b = coeff_ade._b;
          const auto e_cur = ex(i, j, k);
          ex(i, j, k) = cexe(i, j, k) * ex(i, j, k) +
                        cexhy(i, j, k) * (hy(i, j, k) - hy(i, j, k - 1)) +
                        cexhz(i, j, k) * (hz(i, j, k) - hz(i, j - 1, k)) -
                        b * j_sum;

          update_j(i, j, k, coeff_ade, dt, ex(i, j, k), e_cur, j_record._jx);
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

        if (m_index == -1 || _drude_map.find(m_index) == _drude_map.end()) {
          ey(i, j, k) = ceye(i, j, k) * ey(i, j, k) +
                        ceyhz(i, j, k) * (hz(i, j, k) - hz(i - 1, j, k)) +
                        ceyhx(i, j, k) * (hx(i, j, k) - hx(i, j, k - 1));
        } else {
          auto l_index = _drude_map.find(m_index)->second;
          const auto& coeff_ade = _coeff[l_index];
          auto& j_record = _j_record[l_index];

          double j_sum = get_j_sum(i, j, k, coeff_ade, j_record._jy);

          const auto& b = coeff_ade._b;
          const auto e_cur = ey(i, j, k);

          ey(i, j, k) = ceye(i, j, k) * ey(i, j, k) +
                        ceyhz(i, j, k) * (hz(i, j, k) - hz(i - 1, j, k)) +
                        ceyhx(i, j, k) * (hx(i, j, k) - hx(i, j, k - 1)) -
                        b * j_sum;

          update_j(i, j, k, coeff_ade, dt, ey(i, j, k), e_cur, j_record._jy);
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

        if (m_index == -1 || _drude_map.find(m_index) == _drude_map.end()) {
          ez(i, j, k) = ceze(i, j, k) * ez(i, j, k) +
                        cezhx(i, j, k) * (hx(i, j, k) - hx(i, j - 1, k)) +
                        cezhy(i, j, k) * (hy(i, j, k) - hy(i - 1, j, k));
        } else {
          auto l_index = _drude_map.find(m_index)->second;
          const auto& coeff_ade = _coeff[l_index];
          auto& j_record = _j_record[l_index];

          double j_sum = get_j_sum(i, j, k, coeff_ade, j_record._jz);

          const auto& b = coeff_ade._b;
          const auto e_prev_temp = ez(i, j, k);

          ez(i, j, k) = ceze(i, j, k) * ez(i, j, k) +
                        cezhx(i, j, k) * (hx(i, j, k) - hx(i, j - 1, k)) +
                        cezhy(i, j, k) * (hy(i, j, k) - hy(i - 1, j, k)) -
                        b * j_sum;

          update_j(i, j, k, coeff_ade, dt, ez(i, j, k), e_prev_temp,
                   j_record._jz);
        }
      }
    }
  }

  updateEEdge();
}

auto DrudeADEUpdator::updateEEdge() -> void {
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

  const auto& dt = _calculation_param->timeParam()->dt();

  bool contain_xn_edge = containXNEdge();
  bool contain_yn_edge = containYNEdge();
  bool contain_zn_edge = containZNEdge();

  auto get_j_sum = [](std::size_t i, std::size_t j, std::size_t k,
                      const auto& coeff, const auto& j_arr) {
    double j_sum{0};
    auto num_p = coeff._k.size();
    for (decltype(num_p) p = 0; p < num_p; ++p) {
      const auto& k = coeff._k(p);
      const auto& jj = j_arr(i, j, k, p);
      j_sum += ((1 + k)) * jj;
    }
    return 0.5 * j_sum;
  };

  auto update_j = [](std::size_t i, std::size_t j, std::size_t k,
                     const auto& coeff, const auto& dt, const auto& e_next,
                     const auto& e_cur, auto&& jj) {
    auto num_p = coeff._k.size();
    for (decltype(num_p) p = 0; p < num_p; ++p) {
      const auto& k = coeff._k(p);
      const auto& beta = coeff._beta(p);

      jj(i, j, k, p) = k * jj(i, j, k, p) + beta * (e_next + e_cur);
    }
  };

  if (!contain_yn_edge && !contain_zn_edge) {
    auto j = js;
    auto k = ks;
    for (std::size_t i{is}; i < ie; ++i) {
      auto m_index = _grid_space->gridWithMaterial()(i, j, k)->materialIndex();
      if (m_index == -1 || _drude_map.find(m_index) == _drude_map.end()) {
        ex(i, j, k) = cexe(i, j, k) * ex(i, j, k) +
                      cexhy(i, j, k) * (hy(i, j, k) - hy(i, j, k - 1)) +
                      cexhz(i, j, k) * (hz(i, j, k) - hz(i, j - 1, k));
      } else {
        auto l_index = _drude_map.find(m_index)->second;
        const auto& coeff_ade = _coeff[l_index];
        auto& j_record = _j_record[l_index];

        double j_sum = get_j_sum(i, j, k, coeff_ade, j_record._jx);

        const auto& b = coeff_ade._b;
        const auto e_cur = ex(i, j, k);
        ex(i, j, k) = cexe(i, j, k) * ex(i, j, k) +
                      cexhy(i, j, k) * (hy(i, j, k) - hy(i, j, k - 1)) +
                      cexhz(i, j, k) * (hz(i, j, k) - hz(i, j - 1, k)) -
                      b * j_sum;

        update_j(i, j, k, coeff_ade, dt, ex(i, j, k), e_cur, j_record._jx);
      }
    }
  }

  if (!contain_xn_edge && !contain_zn_edge) {
    auto i = is;
    auto k = ks;
    for (std::size_t j{js}; j < je; ++j) {
      auto m_index = _grid_space->gridWithMaterial()(i, j, k)->materialIndex();

      if (m_index == -1 || _drude_map.find(m_index) == _drude_map.end()) {
        ey(i, j, k) = ceye(i, j, k) * ey(i, j, k) +
                      ceyhz(i, j, k) * (hz(i, j, k) - hz(i - 1, j, k)) +
                      ceyhx(i, j, k) * (hx(i, j, k) - hx(i, j, k - 1));
      } else {
        auto l_index = _drude_map.find(m_index)->second;
        const auto& coeff_ade = _coeff[l_index];
        auto& j_record = _j_record[l_index];

        double j_sum = get_j_sum(i, j, k, coeff_ade, j_record._jy);

        const auto& b = coeff_ade._b;
        const auto e_cur = ey(i, j, k);

        ey(i, j, k) = ceye(i, j, k) * ey(i, j, k) +
                      ceyhz(i, j, k) * (hz(i, j, k) - hz(i - 1, j, k)) +
                      ceyhx(i, j, k) * (hx(i, j, k) - hx(i, j, k - 1)) -
                      b * j_sum;

        update_j(i, j, k, coeff_ade, dt, ey(i, j, k), e_cur, j_record._jy);
      }
    }
  }

  if (!contain_xn_edge && !contain_yn_edge) {
    auto i = is;
    auto j = js;
    for (std::size_t k{ks}; k < ke; ++k) {
      auto m_index = _grid_space->gridWithMaterial()(i, j, k)->materialIndex();

      if (m_index == -1 || _drude_map.find(m_index) == _drude_map.end()) {
        ez(i, j, k) = ceze(i, j, k) * ez(i, j, k) +
                      cezhx(i, j, k) * (hx(i, j, k) - hx(i, j - 1, k)) +
                      cezhy(i, j, k) * (hy(i, j, k) - hy(i - 1, j, k));
      } else {
        auto l_index = _drude_map.find(m_index)->second;
        const auto& coeff_ade = _coeff[l_index];
        auto& j_record = _j_record[l_index];

        double j_sum = get_j_sum(i, j, k, coeff_ade, j_record._jz);

        const auto& b = coeff_ade._b;
        const auto e_prev_temp = ez(i, j, k);

        ez(i, j, k) = ceze(i, j, k) * ez(i, j, k) +
                      cezhx(i, j, k) * (hx(i, j, k) - hx(i, j - 1, k)) +
                      cezhy(i, j, k) * (hy(i, j, k) - hy(i - 1, j, k)) -
                      b * j_sum;

        update_j(i, j, k, coeff_ade, dt, ez(i, j, k), e_prev_temp,
                 j_record._jz);
      }
    }
  }

  if (!contain_xn_edge) {
    auto i = is;
    for (std::size_t j{js}; j < je; ++j) {
      for (std::size_t k{ks + 1}; k < ke; ++k) {
        auto m_index =
            _grid_space->gridWithMaterial()(i, j, k)->materialIndex();

        if (m_index == -1 || _drude_map.find(m_index) == _drude_map.end()) {
          ey(i, j, k) = ceye(i, j, k) * ey(i, j, k) +
                        ceyhz(i, j, k) * (hz(i, j, k) - hz(i - 1, j, k)) +
                        ceyhx(i, j, k) * (hx(i, j, k) - hx(i, j, k - 1));
        } else {
          auto l_index = _drude_map.find(m_index)->second;
          const auto& coeff_ade = _coeff[l_index];
          auto& j_record = _j_record[l_index];

          double j_sum = get_j_sum(i, j, k, coeff_ade, j_record._jy);

          const auto& b = coeff_ade._b;
          const auto e_cur = ey(i, j, k);

          ey(i, j, k) = ceye(i, j, k) * ey(i, j, k) +
                        ceyhz(i, j, k) * (hz(i, j, k) - hz(i - 1, j, k)) +
                        ceyhx(i, j, k) * (hx(i, j, k) - hx(i, j, k - 1)) -
                        b * j_sum;

          update_j(i, j, k, coeff_ade, dt, ey(i, j, k), e_cur, j_record._jy);
        }
      }
    }
    for (std::size_t j{js + 1}; j < je; ++j) {
      for (std::size_t k{ks}; k < ke; ++k) {
        auto m_index =
            _grid_space->gridWithMaterial()(i, j, k)->materialIndex();

        if (m_index == -1 || _drude_map.find(m_index) == _drude_map.end()) {
          ez(i, j, k) = ceze(i, j, k) * ez(i, j, k) +
                        cezhx(i, j, k) * (hx(i, j, k) - hx(i, j - 1, k)) +
                        cezhy(i, j, k) * (hy(i, j, k) - hy(i - 1, j, k));
        } else {
          auto l_index = _drude_map.find(m_index)->second;
          const auto& coeff_ade = _coeff[l_index];
          auto& j_record = _j_record[l_index];

          double j_sum = get_j_sum(i, j, k, coeff_ade, j_record._jz);

          const auto& b = coeff_ade._b;
          const auto e_prev_temp = ez(i, j, k);

          ez(i, j, k) = ceze(i, j, k) * ez(i, j, k) +
                        cezhx(i, j, k) * (hx(i, j, k) - hx(i, j - 1, k)) +
                        cezhy(i, j, k) * (hy(i, j, k) - hy(i - 1, j, k)) -
                        b * j_sum;

          update_j(i, j, k, coeff_ade, dt, ez(i, j, k), e_prev_temp,
                   j_record._jz);
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

        if (m_index == -1 || _drude_map.find(m_index) == _drude_map.end()) {
          ez(i, j, k) = ceze(i, j, k) * ez(i, j, k) +
                        cezhx(i, j, k) * (hx(i, j, k) - hx(i, j - 1, k)) +
                        cezhy(i, j, k) * (hy(i, j, k) - hy(i - 1, j, k));
        } else {
          auto l_index = _drude_map.find(m_index)->second;
          const auto& coeff_ade = _coeff[l_index];
          auto& j_record = _j_record[l_index];

          double j_sum = get_j_sum(i, j, k, coeff_ade, j_record._jz);

          const auto& b = coeff_ade._b;
          const auto e_prev_temp = ez(i, j, k);

          ez(i, j, k) = ceze(i, j, k) * ez(i, j, k) +
                        cezhx(i, j, k) * (hx(i, j, k) - hx(i, j - 1, k)) +
                        cezhy(i, j, k) * (hy(i, j, k) - hy(i - 1, j, k)) -
                        b * j_sum;

          update_j(i, j, k, coeff_ade, dt, ez(i, j, k), e_prev_temp,
                   j_record._jz);
        }
      }
    }
    for (std::size_t i{is}; i < ie; ++i) {
      for (std::size_t k{ks + 1}; k < ke; ++k) {
        auto m_index =
            _grid_space->gridWithMaterial()(i, j, k)->materialIndex();
        if (m_index == -1 || _drude_map.find(m_index) == _drude_map.end()) {
          ex(i, j, k) = cexe(i, j, k) * ex(i, j, k) +
                        cexhy(i, j, k) * (hy(i, j, k) - hy(i, j, k - 1)) +
                        cexhz(i, j, k) * (hz(i, j, k) - hz(i, j - 1, k));
        } else {
          auto l_index = _drude_map.find(m_index)->second;
          const auto& coeff_ade = _coeff[l_index];
          auto& j_record = _j_record[l_index];

          double j_sum = get_j_sum(i, j, k, coeff_ade, j_record._jx);

          const auto& b = coeff_ade._b;
          const auto e_cur = ex(i, j, k);
          ex(i, j, k) = cexe(i, j, k) * ex(i, j, k) +
                        cexhy(i, j, k) * (hy(i, j, k) - hy(i, j, k - 1)) +
                        cexhz(i, j, k) * (hz(i, j, k) - hz(i, j - 1, k)) -
                        b * j_sum;

          update_j(i, j, k, coeff_ade, dt, ex(i, j, k), e_cur, j_record._jx);
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
        if (m_index == -1 || _drude_map.find(m_index) == _drude_map.end()) {
          ex(i, j, k) = cexe(i, j, k) * ex(i, j, k) +
                        cexhy(i, j, k) * (hy(i, j, k) - hy(i, j, k - 1)) +
                        cexhz(i, j, k) * (hz(i, j, k) - hz(i, j - 1, k));
        } else {
          auto l_index = _drude_map.find(m_index)->second;
          const auto& coeff_ade = _coeff[l_index];
          auto& j_record = _j_record[l_index];

          double j_sum = get_j_sum(i, j, k, coeff_ade, j_record._jx);

          const auto& b = coeff_ade._b;
          const auto e_cur = ex(i, j, k);
          ex(i, j, k) = cexe(i, j, k) * ex(i, j, k) +
                        cexhy(i, j, k) * (hy(i, j, k) - hy(i, j, k - 1)) +
                        cexhz(i, j, k) * (hz(i, j, k) - hz(i, j - 1, k)) -
                        b * j_sum;

          update_j(i, j, k, coeff_ade, dt, ex(i, j, k), e_cur, j_record._jx);
        }
      }
    }
    for (std::size_t i{is + 1}; i < ie; ++i) {
      for (std::size_t j{js}; j < je; ++j) {
        auto m_index =
            _grid_space->gridWithMaterial()(i, j, k)->materialIndex();

        if (m_index == -1 || _drude_map.find(m_index) == _drude_map.end()) {
          ey(i, j, k) = ceye(i, j, k) * ey(i, j, k) +
                        ceyhz(i, j, k) * (hz(i, j, k) - hz(i - 1, j, k)) +
                        ceyhx(i, j, k) * (hx(i, j, k) - hx(i, j, k - 1));
        } else {
          auto l_index = _drude_map.find(m_index)->second;
          const auto& coeff_ade = _coeff[l_index];
          auto& j_record = _j_record[l_index];

          double j_sum = get_j_sum(i, j, k, coeff_ade, j_record._jy);

          const auto& b = coeff_ade._b;
          const auto e_cur = ey(i, j, k);

          ey(i, j, k) = ceye(i, j, k) * ey(i, j, k) +
                        ceyhz(i, j, k) * (hz(i, j, k) - hz(i - 1, j, k)) +
                        ceyhx(i, j, k) * (hx(i, j, k) - hx(i, j, k - 1)) -
                        b * j_sum;

          update_j(i, j, k, coeff_ade, dt, ey(i, j, k), e_cur, j_record._jy);
        }
      }
    }
  }
}

}  // namespace xfdtd
