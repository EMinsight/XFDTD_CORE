#include "updator/dispersive_material_updator.h"

#include <xfdtd/material/dispersive_material.h>
#include <xfdtd/util/fdtd_basic.h>

#include <memory>
#include <utility>

#include "updator/update_scheme.h"

namespace xfdtd {

LinearDispersiveMaterialADEUpdator::LinearDispersiveMaterialADEUpdator(
    std::vector<std::shared_ptr<Material>> material_arr,
    std::shared_ptr<const GridSpace> grid_space,
    std::shared_ptr<const CalculationParam> calculation_param,
    std::shared_ptr<EMF> emf, IndexTask task)
    : BasicUpdator3D(std::move(grid_space), std::move(calculation_param),
                     std::move(emf), task) {
  init(std::move(material_arr));
}

LinearDispersiveMaterialADEUpdator::ADECorrector::ADECorrector(
    Index num_pole, IndexTask task, Real coeff_j,
    std::shared_ptr<const CalculationParam> calculation_param,
    std::shared_ptr<EMF> emf)
    : _jx{xt::zeros<Real>({task.xRange().size(), task.yRange().size() + 1,
                           task.zRange().size() + 1, num_pole})},
      _jy{xt::zeros<Real>({task.xRange().size() + 1, task.yRange().size(),
                           task.zRange().size() + 1, num_pole})},
      _jz{xt::zeros<Real>({task.xRange().size() + 1, task.yRange().size() + 1,
                           task.zRange().size(), num_pole})},
      _num_pole{num_pole},
      _task{task},
      _coeff_j{coeff_j},
      _calculation_param{std::move(calculation_param)},
      _emf{std::move(emf)} {
  if (num_pole < 1) {
    throw XFDTDLinearDispersiveMaterialException{
        "Number of poles must be greater than 0"};
  }

  if (!task.valid()) {
    throw XFDTDLinearDispersiveMaterialException{"Invalid task"};
  }
}

auto LinearDispersiveMaterialADEUpdator::init(
    std::vector<std::shared_ptr<Material>> material_arr) -> void {
  if (material_arr.empty()) {
    throw XFDTDLinearDispersiveMaterialException{"Material array is empty"};
  }

  _map.resize(material_arr.size(), -1);

  for (Index i{0}; i < material_arr.size(); ++i) {
    const auto& m = material_arr[i];
    if (!m->dispersion()) {
      continue;
    }

    auto dispersive_material =
        std::dynamic_pointer_cast<LinearDispersiveMaterial>(m);
    if (dispersive_material == nullptr) {
      continue;
    }

    auto type = dispersive_material->type();

    switch (type) {
      case LinearDispersiveMaterial::Type::DEBYE: {
        const auto& d =
            std::dynamic_pointer_cast<DebyeMedium>(dispersive_material);

        _ade_correctors.emplace_back(std::make_unique<DebyeADECorrector>(
            *d, _calculation_param, _emf, task()));
        _map[i] = _ade_correctors.size() - 1;
        break;
      }
      case LinearDispersiveMaterial::Type::DRUDE: {
        const auto& d =
            std::dynamic_pointer_cast<DrudeMedium>(dispersive_material);

        _ade_correctors.emplace_back(std::make_unique<DrudeADECorrector>(
            *d, _calculation_param, _emf, task()));
        _map[i] = _ade_correctors.size() - 1;
        break;
      }
      case LinearDispersiveMaterial::Type::LORENTZ: {
        const auto& l =
            std::dynamic_pointer_cast<LorentzMedium>(dispersive_material);

        _ade_correctors.emplace_back(std::make_unique<LorentzADECorrector>(
            *l, _calculation_param, _emf, task()));
        _map[i] = _ade_correctors.size() - 1;
        break;
      }
      default:
        throw XFDTDLinearDispersiveMaterialException{"Invalid material type"};
    }
  }
}

auto LinearDispersiveMaterialADEUpdator::updateE() -> void {
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

        if (m_index == -1 || _map[m_index] == -1) {
          ex(i, j, k) = eNext(cexe(i, j, k), ex(i, j, k), cexhy(i, j, k),
                              hy(i, j, k), hy(i, j, k - 1), cexhz(i, j, k),
                              hz(i, j, k), hz(i, j - 1, k));
          continue;
        }

        _ade_correctors[_map[m_index]]->updateEx(i, j, k);
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

        if (m_index == -1 || _map[m_index] == -1) {
          ey(i, j, k) = eNext(ceye(i, j, k), ey(i, j, k), ceyhz(i, j, k),
                              hz(i, j, k), hz(i - 1, j, k), ceyhx(i, j, k),
                              hx(i, j, k), hx(i, j, k - 1));
          continue;
        }

        _ade_correctors[_map[m_index]]->updateEy(i, j, k);
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

        if (m_index == -1 || _map[m_index] == -1) {
          ez(i, j, k) = eNext(ceze(i, j, k), ez(i, j, k), cezhx(i, j, k),
                              hx(i, j, k), hx(i, j - 1, k), cezhy(i, j, k),
                              hy(i, j, k), hy(i - 1, j, k));
          continue;
        }

        _ade_correctors[_map[m_index]]->updateEz(i, j, k);
      }
    }
  }

  updateEEdge();
}

auto LinearDispersiveMaterialADEUpdator::updateEEdge() -> void {
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

  bool contain_xn_edge = containXNEdge();
  bool contain_yn_edge = containYNEdge();
  bool contain_zn_edge = containZNEdge();

  if (!contain_yn_edge && !contain_zn_edge) {
    auto j = js;
    auto k = ks;
    for (std::size_t i{is}; i < ie; ++i) {
      auto m_index = _grid_space->gridWithMaterial()(i, j, k)->materialIndex();
      if (m_index == -1 || _map[m_index] == -1) {
        ex(i, j, k) = eNext(cexe(i, j, k), ex(i, j, k), cexhy(i, j, k),
                            hy(i, j, k), hy(i, j, k - 1), cexhz(i, j, k),
                            hz(i, j, k), hz(i, j - 1, k));
        continue;
      }

      _ade_correctors[_map[m_index]]->updateEx(i, j, k);
    }
  }

  if (!contain_xn_edge && !contain_zn_edge) {
    auto i = is;
    auto k = ks;
    for (std::size_t j{js}; j < je; ++j) {
      auto m_index = _grid_space->gridWithMaterial()(i, j, k)->materialIndex();

      if (m_index == -1 || _map[m_index] == -1) {
        ey(i, j, k) = eNext(ceye(i, j, k), ey(i, j, k), ceyhz(i, j, k),
                            hz(i, j, k), hz(i - 1, j, k), ceyhx(i, j, k),
                            hx(i, j, k), hx(i, j, k - 1));
        continue;
      }

      _ade_correctors[_map[m_index]]->updateEy(i, j, k);
    }
  }

  if (!contain_xn_edge && !contain_yn_edge) {
    auto i = is;
    auto j = js;
    for (std::size_t k{ks}; k < ke; ++k) {
      auto m_index = _grid_space->gridWithMaterial()(i, j, k)->materialIndex();

      if (m_index == -1 || _map[m_index] == -1) {
        ez(i, j, k) = eNext(ceze(i, j, k), ez(i, j, k), cezhx(i, j, k),
                            hx(i, j, k), hx(i, j - 1, k), cezhy(i, j, k),
                            hy(i, j, k), hy(i - 1, j, k));
        continue;
      }

      _ade_correctors[_map[m_index]]->updateEz(i, j, k);
    }
  }

  if (!contain_xn_edge) {
    auto i = is;
    for (std::size_t j{js}; j < je; ++j) {
      for (std::size_t k{ks + 1}; k < ke; ++k) {
        auto m_index =
            _grid_space->gridWithMaterial()(i, j, k)->materialIndex();

        if (m_index == -1 || _map[m_index] == -1) {
          ey(i, j, k) = eNext(ceye(i, j, k), ey(i, j, k), ceyhz(i, j, k),
                              hz(i, j, k), hz(i - 1, j, k), ceyhx(i, j, k),
                              hx(i, j, k), hx(i, j, k - 1));
          continue;
        }

        _ade_correctors[_map[m_index]]->updateEy(i, j, k);
      }
    }
    for (std::size_t j{js + 1}; j < je; ++j) {
      for (std::size_t k{ks}; k < ke; ++k) {
        auto m_index =
            _grid_space->gridWithMaterial()(i, j, k)->materialIndex();

        if (m_index == -1 || _map[m_index] == -1) {
          ez(i, j, k) = eNext(ceze(i, j, k), ez(i, j, k), cezhx(i, j, k),
                              hx(i, j, k), hx(i, j - 1, k), cezhy(i, j, k),
                              hy(i, j, k), hy(i - 1, j, k));
          continue;
        }

        _ade_correctors[_map[m_index]]->updateEz(i, j, k);
      }
    }
  }

  if (!contain_yn_edge) {
    auto j = js;
    for (std::size_t i{is + 1}; i < ie; ++i) {
      for (std::size_t k{ks}; k < ke; ++k) {
        auto m_index =
            _grid_space->gridWithMaterial()(i, j, k)->materialIndex();

        if (m_index == -1 || _map[m_index] == -1) {
          ez(i, j, k) = eNext(ceze(i, j, k), ez(i, j, k), cezhx(i, j, k),
                              hx(i, j, k), hx(i, j - 1, k), cezhy(i, j, k),
                              hy(i, j, k), hy(i - 1, j, k));
          continue;
        }

        _ade_correctors[_map[m_index]]->updateEz(i, j, k);
      }
    }
    for (std::size_t i{is}; i < ie; ++i) {
      for (std::size_t k{ks + 1}; k < ke; ++k) {
        auto m_index =
            _grid_space->gridWithMaterial()(i, j, k)->materialIndex();
        if (m_index == -1 || _map[m_index] == -1) {
          ex(i, j, k) = eNext(cexe(i, j, k), ex(i, j, k), cexhy(i, j, k),
                              hy(i, j, k), hy(i, j, k - 1), cexhz(i, j, k),
                              hz(i, j, k), hz(i, j - 1, k));
          continue;
        }

        _ade_correctors[_map[m_index]]->updateEx(i, j, k);
      }
    }
  }

  if (!contain_zn_edge) {
    auto k = ks;
    for (std::size_t i{is}; i < ie; ++i) {
      for (std::size_t j{js + 1}; j < je; ++j) {
        auto m_index =
            _grid_space->gridWithMaterial()(i, j, k)->materialIndex();
        if (m_index == -1 || _map[m_index] == -1) {
          ex(i, j, k) = eNext(cexe(i, j, k), ex(i, j, k), cexhy(i, j, k),
                              hy(i, j, k), hy(i, j, k - 1), cexhz(i, j, k),
                              hz(i, j, k), hz(i, j - 1, k));
          continue;
        }

        _ade_correctors[_map[m_index]]->updateEx(i, j, k);
      }
    }
    for (std::size_t i{is + 1}; i < ie; ++i) {
      for (std::size_t j{js}; j < je; ++j) {
        auto m_index =
            _grid_space->gridWithMaterial()(i, j, k)->materialIndex();

        if (m_index == -1 || _map[m_index] == -1) {
          ey(i, j, k) = eNext(ceye(i, j, k), ey(i, j, k), ceyhz(i, j, k),
                              hz(i, j, k), hz(i - 1, j, k), ceyhx(i, j, k),
                              hx(i, j, k), hx(i, j, k - 1));
          continue;
        }

        _ade_correctors[_map[m_index]]->updateEy(i, j, k);
      }
    }
  }
}

}  // namespace xfdtd
