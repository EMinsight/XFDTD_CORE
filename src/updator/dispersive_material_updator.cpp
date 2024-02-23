#include "updator/dispersive_material_updator.h"

#include <cassert>
#include <cstddef>
#include <memory>
#include <utility>
#include <vector>

#include "updator/basic_updator.h"
#include "xfdtd/grid_space/grid_space.h"
#include "xfdtd/material/dispersive_material.h"

namespace xfdtd {

LinearDispersiveMaterialADEUpdator::LinearDispersiveMaterialADEUpdator(
    std::shared_ptr<const GridSpace> grid_space,
    std::shared_ptr<const CalculationParam> calculation_param,
    std::shared_ptr<EMF> emf)
    : BasicUpdator(std::move(grid_space), std::move(calculation_param),
                   std::move(emf)) {}

LorentzADEUpdator::LorentzADEUpdator(
    std::shared_ptr<const GridSpace> grid_space,
    std::shared_ptr<const CalculationParam> calculation_param,
    std::shared_ptr<EMF> emf)
    : LinearDispersiveMaterialADEUpdator{
          std::move(grid_space), std::move(calculation_param), std::move(emf)} {
  assert(_grid_space->type() != GridSpace::Type::NONUNIFORM);
  init();
}

void LorentzADEUpdator::updateE() {
  const auto nx{_grid_space->sizeX()};
  const auto ny{_grid_space->sizeY()};
  const auto nz{_grid_space->sizeZ()};

  // FIXME: don't use xt::view!! too slow!!
  // is there ant nice way to do this?
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

  for (std::size_t i{0}; i < nx; ++i) {
    for (std::size_t j{1}; j < ny; ++j) {
      for (std::size_t k{1}; k < nz; ++k) {
        auto m_index = _grid_space->grid()(i, j, k)->materialIndex();
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

  for (std::size_t i{1}; i < nx; ++i) {
    for (std::size_t j{0}; j < ny; ++j) {
      for (std::size_t k{1}; k < nz; ++k) {
        auto m_index = _grid_space->grid()(i, j, k)->materialIndex();

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

  for (std::size_t i{1}; i < nx; ++i) {
    for (std::size_t j{1}; j < ny; ++j) {
      for (std::size_t k{0}; k < nz; ++k) {
        auto m_index = _grid_space->grid()(i, j, k)->materialIndex();

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
}

void LorentzADEUpdator::init() { allocate(); }

void LorentzADEUpdator::allocate() {
  const auto& nx = _grid_space->sizeX();
  const auto& ny = _grid_space->sizeY();
  const auto& nz = _grid_space->sizeZ();

  _emf->allocateExPrev(nx, ny + 1, nz + 1);
  _emf->allocateEyPrev(nx + 1, ny, nz + 1);
  _emf->allocateEzPrev(nx + 1, ny + 1, nz);

  const auto& num_material =
      _calculation_param->materialParam()->materialArray().size();
  std::size_t m_index = 0;
  std::size_t temp_i = 0;
  _lorentz_mediums = {};
  for (const auto& m : _calculation_param->materialParam()->materialArray()) {
    if (!m->dispersion()) {
      ++m_index;
      continue;
    }

    auto linear_dispersive_material =
        std::dynamic_pointer_cast<LinearDispersiveMaterial>(m);
    if (linear_dispersive_material == nullptr) {
      ++m_index;
      continue;
    }

    auto lorentz_medium =
        std::dynamic_pointer_cast<LorentzMedium>(linear_dispersive_material);
    if (lorentz_medium == nullptr) {
      ++m_index;
      continue;
    }

    _lorentz_map.insert({m_index, temp_i});
    _lorentz_mediums.emplace_back(lorentz_medium);
    ++m_index;
    ++temp_i;
  }

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

}  // namespace xfdtd
