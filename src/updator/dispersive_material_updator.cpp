#include "updator/dispersive_material_updator.h"

#include <cstddef>
#include <memory>
#include <utility>
#include <vector>

#include "updator/basic_updator.h"

namespace xfdtd {

template <typename T>
static void handleDispersiveMaterialADEUpdator(
    std::unordered_map<std::size_t, std::size_t>& map,
    std::vector<std::shared_ptr<T>>& dispersion_arr,
    const std::vector<std::shared_ptr<Material>>& material_arr) {
  std::size_t m_index = 0;
  std::size_t temp_i = 0;
  for (const auto& m : material_arr) {
    if (!m->dispersion()) {
      ++m_index;
      continue;
    }

    auto dispersive_material = std::dynamic_pointer_cast<T>(m);
    if (dispersive_material == nullptr) {
      ++m_index;
      continue;
    }

    map.insert({m_index, temp_i});
    dispersion_arr.emplace_back(dispersive_material);
    ++m_index;
    ++temp_i;
  }
}

LinearDispersiveMaterialADEUpdator::LinearDispersiveMaterialADEUpdator(
    std::shared_ptr<const GridSpace> grid_space,
    std::shared_ptr<const CalculationParam> calculation_param,
    std::shared_ptr<EMF> emf, Divider::IndexTask task)
    : BasicUpdator(std::move(grid_space), std::move(calculation_param),
                   std::move(emf), std::move(task)) {
  if (_grid_space->type() != GridSpace::Type::UNIFORM) {
    throw std::runtime_error("Grid space type must be uniform");
  }
}

LorentzADEUpdator::LorentzADEUpdator(
    std::shared_ptr<const GridSpace> grid_space,
    std::shared_ptr<const CalculationParam> calculation_param,
    std::shared_ptr<EMF> emf, Divider::IndexTask task)
    : LinearDispersiveMaterialADEUpdator{std::move(grid_space),
                                         std::move(calculation_param),
                                         std::move(emf), std::move(task)} {
  init();
}

void LorentzADEUpdator::updateE() {
  const auto is = task()._x_range[0];
  const auto ie = task()._x_range[1];
  const auto js = task()._y_range[0];
  const auto je = task()._y_range[1];
  const auto ks = task()._z_range[0];
  const auto ke = task()._z_range[1];

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
  auto& ex_prev = _emf->exPrev();
  auto& ey_prev = _emf->eyPrev();
  auto& ez_prev = _emf->ezPrev();

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

  for (std::size_t i{is}; i < ie; ++i) {
    for (std::size_t j{js + 1}; j < je; ++j) {
      for (std::size_t k{ks + 1}; k < ke; ++k) {
        auto m_index = _grid_space->grid()(i, j, k)->materialIndex();
        if (m_index == -1 || _lorentz_map.find(m_index) == _lorentz_map.end()) {
          ex(i, j, k) = cexe(i, j, k) * ex(i, j, k) +
                        cexhy(i, j, k) * (hy(i, j, k) - hy(i, j, k - 1)) +
                        cexhz(i, j, k) * (hz(i, j, k) - hz(i, j - 1, k));
        } else {
          auto l_index = _lorentz_map.find(m_index)->second;
          const auto& coeff_ade = _coeff[l_index];
          auto& j_record = _j_record[l_index];

          auto local_i = i - is;
          auto local_j = j - js;
          auto local_k = k - ks;
          double j_sum = get_j_sum(local_i, local_j, local_k, coeff_ade,
                                   j_record._jx, j_record._jx_prev);

          const auto& c1 = coeff_ade._c1;
          const auto& c2 = coeff_ade._c2;
          const auto& c3 = coeff_ade._c3;
          const auto e_prev_temp = ex(i, j, k);
          ex(i, j, k) = c1 * ex_prev(i, j, k) + cexe(i, j, k) * ex(i, j, k) +
                        cexhy(i, j, k) * (hy(i, j, k) - hy(i, j, k - 1)) +
                        cexhz(i, j, k) * (hz(i, j, k) - hz(i, j - 1, k)) -
                        c3 * j_sum;

          update_j(local_i, local_j, local_k, coeff_ade, dt, ex(i, j, k),
                   ex_prev(i, j, k), j_record._jx, j_record._jx_prev);
          ex_prev(i, j, k) = e_prev_temp;
        }
      }
    }
  }

  for (std::size_t i{is + 1}; i < ie; ++i) {
    for (std::size_t j{js}; j < je; ++j) {
      for (std::size_t k{ks + 1}; k < ke; ++k) {
        auto m_index = _grid_space->grid()(i, j, k)->materialIndex();

        if (m_index == -1 || _lorentz_map.find(m_index) == _lorentz_map.end()) {
          ey(i, j, k) = ceye(i, j, k) * ey(i, j, k) +
                        ceyhz(i, j, k) * (hz(i, j, k) - hz(i - 1, j, k)) +
                        ceyhx(i, j, k) * (hx(i, j, k) - hx(i, j, k - 1));
        } else {
          auto l_index = _lorentz_map.find(m_index)->second;
          const auto& coeff_ade = _coeff[l_index];
          auto& j_record = _j_record[l_index];

          auto local_i = i - is;
          auto local_j = j - js;
          auto local_k = k - ks;
          double j_sum = get_j_sum(local_i, local_j, local_k, coeff_ade,
                                   j_record._jy, j_record._jy_prev);

          const auto& c1 = coeff_ade._c1;
          const auto& c2 = coeff_ade._c2;
          const auto& c3 = coeff_ade._c3;
          const auto e_prev_temp = ey(i, j, k);

          ey(i, j, k) = c1 * ey_prev(i, j, k) + ceye(i, j, k) * ey(i, j, k) +
                        ceyhz(i, j, k) * (hz(i, j, k) - hz(i - 1, j, k)) +
                        ceyhx(i, j, k) * (hx(i, j, k) - hx(i, j, k - 1)) -
                        c3 * j_sum;

          update_j(local_i, local_j, local_k, coeff_ade, dt, ey(i, j, k),
                   ey_prev(i, j, k), j_record._jy, j_record._jy_prev);
          ey_prev(i, j, k) = e_prev_temp;
        }
      }
    }
  }

  for (std::size_t i{is + 1}; i < ie; ++i) {
    for (std::size_t j{js + 1}; j < je; ++j) {
      for (std::size_t k{ks}; k < ke; ++k) {
        auto m_index = _grid_space->grid()(i, j, k)->materialIndex();

        if (m_index == -1 || _lorentz_map.find(m_index) == _lorentz_map.end()) {
          ez(i, j, k) = ceze(i, j, k) * ez(i, j, k) +
                        cezhx(i, j, k) * (hx(i, j, k) - hx(i, j - 1, k)) +
                        cezhy(i, j, k) * (hy(i, j, k) - hy(i - 1, j, k));
        } else {
          auto l_index = _lorentz_map.find(m_index)->second;
          const auto& coeff_ade = _coeff[l_index];
          auto& j_record = _j_record[l_index];

          auto local_i = i - is;
          auto local_j = j - js;
          auto local_k = k - ks;
          double j_sum = get_j_sum(local_i, local_j, local_k, coeff_ade,
                                   j_record._jz, j_record._jz_prev);

          const auto& c1 = coeff_ade._c1;
          const auto& c2 = coeff_ade._c2;
          const auto& c3 = coeff_ade._c3;
          const auto e_prev_temp = ez(i, j, k);

          ez(i, j, k) = c1 * ez_prev(i, j, k) + ceze(i, j, k) * ez(i, j, k) +
                        cezhx(i, j, k) * (hx(i, j, k) - hx(i, j - 1, k)) +
                        cezhy(i, j, k) * (hy(i, j, k) - hy(i - 1, j, k)) -
                        c3 * j_sum;

          update_j(local_i, local_j, local_k, coeff_ade, dt, ez(i, j, k),
                   ez_prev(i, j, k), j_record._jz, j_record._jz_prev);
          ez_prev(i, j, k) = e_prev_temp;
        }
      }
    }
  }
}

void LorentzADEUpdator::init() { allocate(); }

void LorentzADEUpdator::allocate() {
  _lorentz_map = {};
  _lorentz_mediums = {};
  handleDispersiveMaterialADEUpdator(
      _lorentz_map, _lorentz_mediums,
      _calculation_param->materialParam()->materialArray());

  const auto& nx = task()._x_range[1] - task()._x_range[0];
  const auto& ny = task()._y_range[1] - task()._y_range[0];
  const auto& nz = task()._z_range[1] - task()._z_range[0];
  const auto& nm = _lorentz_mediums.size();
  _coeff.resize(nm);
  _j_record.resize(nm);
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

DrudeADEUpdator::DrudeADEUpdator(
    std::shared_ptr<const GridSpace> grid_space,
    std::shared_ptr<const CalculationParam> calculation_param,
    std::shared_ptr<EMF> emf, Divider::IndexTask task)
    : LinearDispersiveMaterialADEUpdator{std::move(grid_space),
                                         std::move(calculation_param),
                                         std::move(emf), std::move(task)} {
  init();
}

void DrudeADEUpdator::init() { allocate(); }

void DrudeADEUpdator::allocate() {
  const auto& nx = task()._x_range[1] - task()._x_range[0];
  const auto& ny = task()._y_range[1] - task()._y_range[0];
  const auto& nz = task()._z_range[1] - task()._z_range[0];

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
  const auto is = task()._x_range[0];
  const auto ie = task()._x_range[1];
  const auto js = task()._y_range[0];
  const auto je = task()._y_range[1];
  const auto ks = task()._z_range[0];
  const auto ke = task()._z_range[1];

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

  for (std::size_t i{is}; i < ie; ++i) {
    for (std::size_t j{js + 1}; j < je; ++j) {
      for (std::size_t k{ks + 1}; k < ke; ++k) {
        auto m_index = _grid_space->grid()(i, j, k)->materialIndex();
        if (m_index == -1 || _drude_map.find(m_index) == _drude_map.end()) {
          ex(i, j, k) = cexe(i, j, k) * ex(i, j, k) +
                        cexhy(i, j, k) * (hy(i, j, k) - hy(i, j, k - 1)) +
                        cexhz(i, j, k) * (hz(i, j, k) - hz(i, j - 1, k));
        } else {
          auto l_index = _drude_map.find(m_index)->second;
          const auto& coeff_ade = _coeff[l_index];
          auto& j_record = _j_record[l_index];

          auto local_i = i - is;
          auto local_j = j - js;
          auto local_k = k - ks;
          double j_sum =
              get_j_sum(local_i, local_j, local_k, coeff_ade, j_record._jx);

          const auto& b = coeff_ade._b;
          const auto e_cur = ex(i, j, k);
          ex(i, j, k) = cexe(i, j, k) * ex(i, j, k) +
                        cexhy(i, j, k) * (hy(i, j, k) - hy(i, j, k - 1)) +
                        cexhz(i, j, k) * (hz(i, j, k) - hz(i, j - 1, k)) -
                        b * j_sum;

          update_j(local_i, local_i, local_k, coeff_ade, dt, ex(i, j, k), e_cur,
                   j_record._jx);
        }
      }
    }
  }

  for (std::size_t i{is + 1}; i < ie; ++i) {
    for (std::size_t j{js}; j < je; ++j) {
      for (std::size_t k{ks + 1}; k < ke; ++k) {
        auto m_index = _grid_space->grid()(i, j, k)->materialIndex();

        if (m_index == -1 || _drude_map.find(m_index) == _drude_map.end()) {
          ey(i, j, k) = ceye(i, j, k) * ey(i, j, k) +
                        ceyhz(i, j, k) * (hz(i, j, k) - hz(i - 1, j, k)) +
                        ceyhx(i, j, k) * (hx(i, j, k) - hx(i, j, k - 1));
        } else {
          auto l_index = _drude_map.find(m_index)->second;
          const auto& coeff_ade = _coeff[l_index];
          auto& j_record = _j_record[l_index];

          auto local_i = i - is;
          auto local_j = j - js;
          auto local_k = k - ks;
          double j_sum =
              get_j_sum(local_i, local_j, local_k, coeff_ade, j_record._jy);

          const auto& b = coeff_ade._b;
          const auto e_cur = ey(i, j, k);

          ey(i, j, k) = ceye(i, j, k) * ey(i, j, k) +
                        ceyhz(i, j, k) * (hz(i, j, k) - hz(i - 1, j, k)) +
                        ceyhx(i, j, k) * (hx(i, j, k) - hx(i, j, k - 1)) -
                        b * j_sum;

          update_j(local_i, local_j, local_k, coeff_ade, dt, ey(i, j, k), e_cur,
                   j_record._jy);
        }
      }
    }
  }

  for (std::size_t i{is + 1}; i < ie; ++i) {
    for (std::size_t j{js + 1}; j < je; ++j) {
      for (std::size_t k{ks}; k < ke; ++k) {
        auto m_index = _grid_space->grid()(i, j, k)->materialIndex();

        if (m_index == -1 || _drude_map.find(m_index) == _drude_map.end()) {
          ez(i, j, k) = ceze(i, j, k) * ez(i, j, k) +
                        cezhx(i, j, k) * (hx(i, j, k) - hx(i, j - 1, k)) +
                        cezhy(i, j, k) * (hy(i, j, k) - hy(i - 1, j, k));
        } else {
          auto l_index = _drude_map.find(m_index)->second;
          const auto& coeff_ade = _coeff[l_index];
          auto& j_record = _j_record[l_index];

          auto local_i = i - is;
          auto local_j = j - js;
          auto local_k = k - ks;
          double j_sum =
              get_j_sum(local_i, local_j, local_k, coeff_ade, j_record._jz);

          const auto& b = coeff_ade._b;
          const auto e_prev_temp = ez(i, j, k);

          ez(i, j, k) = ceze(i, j, k) * ez(i, j, k) +
                        cezhx(i, j, k) * (hx(i, j, k) - hx(i, j - 1, k)) +
                        cezhy(i, j, k) * (hy(i, j, k) - hy(i - 1, j, k)) -
                        b * j_sum;

          update_j(local_i, local_j, local_k, coeff_ade, dt, ez(i, j, k),
                   e_prev_temp, j_record._jz);
        }
      }
    }
  }
}

DebyeADEUpdator::DebyeADEUpdator(
    std::shared_ptr<const GridSpace> grid_space,
    std::shared_ptr<const CalculationParam> calculation_param,
    std::shared_ptr<EMF> emf, Divider::IndexTask task)
    : LinearDispersiveMaterialADEUpdator{std::move(grid_space),
                                         std::move(calculation_param),
                                         std::move(emf), std::move(task)} {
  init();
}

void DebyeADEUpdator::init() {
  const auto& nx = task()._x_range[1] - task()._x_range[0];
  const auto& ny = task()._y_range[1] - task()._y_range[0];
  const auto& nz = task()._z_range[1] - task()._z_range[0];

  _debye_map = {};
  _debye_mediums = {};
  handleDispersiveMaterialADEUpdator(
      _debye_map, _debye_mediums,
      _calculation_param->materialParam()->materialArray());

  const auto& nm = _debye_mediums.size();
  _coeff.resize(nm);
  _j_record.resize(nm);
  for (std::size_t i{0}; i < nm; ++i) {
    _coeff[i] = _debye_mediums[i]->coeffForADE();
    auto num_p = _debye_mediums[i]->numberOfPoles();
    _j_record[i]._jx = xt::zeros<double>({nx, ny + 1, nz + 1, num_p});
    _j_record[i]._jy = xt::zeros<double>({nx + 1, ny, nz + 1, num_p});
    _j_record[i]._jz = xt::zeros<double>({nx + 1, ny + 1, nz, num_p});
  }
}

void DebyeADEUpdator::updateE() {
  const auto is = task()._x_range[0];
  const auto ie = task()._x_range[1];
  const auto js = task()._y_range[0];
  const auto je = task()._y_range[1];
  const auto ks = task()._z_range[0];
  const auto ke = task()._z_range[1];
  const auto dt = _calculation_param->timeParam()->dt();

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

      jj(i, j, k, p) = k * jj(i, j, k, p) + beta * (e_next + e_cur) / dt;
    }
  };

  for (std::size_t i{is}; i < ie; ++i) {
    for (std::size_t j{js + 1}; j < je; ++j) {
      for (std::size_t k{ks + 1}; k < ke; ++k) {
        auto m_index = _grid_space->grid()(i, j, k)->materialIndex();
        if (m_index == -1 || _debye_map.find(m_index) == _debye_map.end()) {
          ex(i, j, k) = cexe(i, j, k) * ex(i, j, k) +
                        cexhy(i, j, k) * (hy(i, j, k) - hy(i, j, k - 1)) +
                        cexhz(i, j, k) * (hz(i, j, k) - hz(i, j - 1, k));
        } else {
          auto l_index = _debye_map.find(m_index)->second;
          const auto& coeff_ade = _coeff[l_index];
          auto& j_record = _j_record[l_index];

          auto local_i = i - is;
          auto local_j = j - js;
          auto local_k = k - ks;
          double j_sum =
              get_j_sum(local_i, local_j, local_k, coeff_ade, j_record._jx);

          const auto& b = coeff_ade._b;
          const auto e_cur = ex(i, j, k);
          ex(i, j, k) = cexe(i, j, k) * ex(i, j, k) +
                        cexhy(i, j, k) * (hy(i, j, k) - hy(i, j, k - 1)) +
                        cexhz(i, j, k) * (hz(i, j, k) - hz(i, j - 1, k)) -
                        b * j_sum;

          update_j(local_i, local_j, local_k, coeff_ade, dt, ex(i, j, k), e_cur,
                   j_record._jx);
        }
      }
    }
  }

  for (std::size_t i{is + 1}; i < ie; ++i) {
    for (std::size_t j{js}; j < je; ++j) {
      for (std::size_t k{ks + 1}; k < ke; ++k) {
        auto m_index = _grid_space->grid()(i, j, k)->materialIndex();

        if (m_index == -1 || _debye_map.find(m_index) == _debye_map.end()) {
          ey(i, j, k) = ceye(i, j, k) * ey(i, j, k) +
                        ceyhz(i, j, k) * (hz(i, j, k) - hz(i - 1, j, k)) +
                        ceyhx(i, j, k) * (hx(i, j, k) - hx(i, j, k - 1));
        } else {
          auto l_index = _debye_map.find(m_index)->second;
          const auto& coeff_ade = _coeff[l_index];
          auto& j_record = _j_record[l_index];

          auto local_i = i - is;
          auto local_j = j - js;
          auto local_k = k - ks;
          double j_sum =
              get_j_sum(local_i, local_j, local_k, coeff_ade, j_record._jy);

          const auto& b = coeff_ade._b;
          const auto e_cur = ey(i, j, k);

          ey(i, j, k) = ceye(i, j, k) * ey(i, j, k) +
                        ceyhz(i, j, k) * (hz(i, j, k) - hz(i - 1, j, k)) +
                        ceyhx(i, j, k) * (hx(i, j, k) - hx(i, j, k - 1)) -
                        b * j_sum;

          update_j(local_i, local_j, local_k, coeff_ade, dt, ey(i, j, k), e_cur,
                   j_record._jy);
        }
      }
    }
  }

  for (std::size_t i{is + 1}; i < ie; ++i) {
    for (std::size_t j{js + 1}; j < je; ++j) {
      for (std::size_t k{ks}; k < ke; ++k) {
        auto m_index = _grid_space->grid()(i, j, k)->materialIndex();

        if (m_index == -1 || _debye_map.find(m_index) == _debye_map.end()) {
          ez(i, j, k) = ceze(i, j, k) * ez(i, j, k) +
                        cezhx(i, j, k) * (hx(i, j, k) - hx(i, j - 1, k)) +
                        cezhy(i, j, k) * (hy(i, j, k) - hy(i - 1, j, k));
        } else {
          auto l_index = _debye_map.find(m_index)->second;
          const auto& coeff_ade = _coeff[l_index];
          auto& j_record = _j_record[l_index];

          auto local_i = i - is;
          auto local_j = j - js;
          auto local_k = k - ks;
          double j_sum =
              get_j_sum(local_i, local_j, local_k, coeff_ade, j_record._jz);

          const auto& b = coeff_ade._b;
          const auto e_prev_temp = ez(i, j, k);

          ez(i, j, k) = ceze(i, j, k) * ez(i, j, k) +
                        cezhx(i, j, k) * (hx(i, j, k) - hx(i, j - 1, k)) +
                        cezhy(i, j, k) * (hy(i, j, k) - hy(i - 1, j, k)) -
                        b * j_sum;

          update_j(local_i, local_j, local_k, coeff_ade, dt, ez(i, j, k),
                   e_prev_temp, j_record._jz);
        }
      }
    }
  }
}

}  // namespace xfdtd
