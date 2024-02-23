#include "updator/basic_updator.h"

#include <xtensor.hpp>

namespace xfdtd {

BasicUpdator::BasicUpdator(
    std::shared_ptr<const GridSpace> grid_space,
    std::shared_ptr<const CalculationParam> calculation_param,
    std::shared_ptr<EMF> emf)
    : Updator(std::move(grid_space), std::move(calculation_param),
              std::move(emf)) {}

void BasicUpdator::updateH() {
  auto nx{_grid_space->sizeX()};
  auto ny{_grid_space->sizeY()};
  auto nz{_grid_space->sizeZ()};

  auto& chxh{_calculation_param->fdtdCoefficient()->chxh()};
  auto& chyh{_calculation_param->fdtdCoefficient()->chyh()};
  auto& chzh{_calculation_param->fdtdCoefficient()->chzh()};
  auto& chxey{_calculation_param->fdtdCoefficient()->chxey()};
  auto& chxez{_calculation_param->fdtdCoefficient()->chxez()};
  auto& chyex{_calculation_param->fdtdCoefficient()->chyex()};
  auto& chyez{_calculation_param->fdtdCoefficient()->chyez()};
  auto& chzex{_calculation_param->fdtdCoefficient()->chzex()};
  auto& chzey{_calculation_param->fdtdCoefficient()->chzey()};

  auto& hx{_emf->hx()};
  auto& hy{_emf->hy()};
  auto& hz{_emf->hz()};
  auto& ex{_emf->ex()};
  auto& ey{_emf->ey()};
  auto& ez{_emf->ez()};

  for (std::size_t i{0}; i < nx; ++i) {
    for (std::size_t j{0}; j < ny; ++j) {
      for (std::size_t k{0}; k < nz; ++k) {
        hx(i, j, k) = chxh(i, j, k) * hx(i, j, k) +
                      chxez(i, j, k) * (ez(i, j + 1, k) - ez(i, j, k)) +
                      chxey(i, j, k) * (ey(i, j, k + 1) - ey(i, j, k));

        hy(i, j, k) = chyh(i, j, k) * hy(i, j, k) +
                      chyez(i, j, k) * (ez(i + 1, j, k) - ez(i, j, k)) +
                      chyex(i, j, k) * (ex(i, j, k + 1) - ex(i, j, k));

        hz(i, j, k) = chzh(i, j, k) * hz(i, j, k) +
                      chzex(i, j, k) * (ex(i, j + 1, k) - ex(i, j, k)) +
                      chzey(i, j, k) * (ey(i + 1, j, k) - ey(i, j, k));
      }
    }
  }
}

BasicUpdatorTEM::BasicUpdatorTEM(
    std::shared_ptr<const GridSpace> grid_space,
    std::shared_ptr<const CalculationParam> calculation_param,
    std::shared_ptr<EMF> emf)
    : BasicUpdator(std::move(grid_space), std::move(calculation_param),
                   std::move(emf)) {}

void BasicUpdatorTEM::updateE() {
  auto nz{_grid_space->hNodeZ().size()};

  auto& cexe{_calculation_param->fdtdCoefficient()->cexe()};
  auto& cexhy{_calculation_param->fdtdCoefficient()->cexhy()};

  auto& hy{_emf->hy()};
  auto& ex{_emf->ex()};

  for (std::size_t k{0}; k < nz; ++k) {
    ex(0, 0, k) = cexe(0, 0, k) * ex(0, 0, k) +
                  cexhy(0, 0, k) * (hy(0, 0, k) - hy(0, 0, k - 1));
  }
}

BasicUpdatorTE::BasicUpdatorTE(
    std::shared_ptr<const GridSpace> grid_space,
    std::shared_ptr<const CalculationParam> calculation_param,
    std::shared_ptr<EMF> emf)
    : BasicUpdator(std::move(grid_space), std::move(calculation_param),
                   std::move(emf)) {}

void BasicUpdatorTE::updateE() {
  auto nx{_grid_space->hNodeX().size()};
  auto ny{_grid_space->hNodeY().size()};
  auto nz{_grid_space->hNodeZ().size()};

  auto& ceze{_calculation_param->fdtdCoefficient()->ceze()};
  auto& cezhx{_calculation_param->fdtdCoefficient()->cezhx()};
  auto& cezhy{_calculation_param->fdtdCoefficient()->cezhy()};

  auto& hx{_emf->hx()};
  auto& hy{_emf->hy()};
  auto& ez{_emf->ez()};

  assert(nz == 1);

  for (std::size_t i{1}; i < nx; ++i) {
    for (std::size_t j{1}; j < ny; ++j) {
      for (std::size_t k{0}; k < nz; ++k) {
        ez(i, j, k) = ceze(i, j, k) * ez(i, j, k) +
                      cezhx(i, j, k) * (hx(i, j, k) - hx(i, j - 1, k)) +
                      cezhy(i, j, k) * (hy(i, j, k) - hy(i - 1, j, k));
      }
    }
  }
}

BasicUpdator3D::BasicUpdator3D(
    std::shared_ptr<const GridSpace> grid_space,
    std::shared_ptr<const CalculationParam> calculation_param,
    std::shared_ptr<EMF> emf)
    : BasicUpdator(std::move(grid_space), std::move(calculation_param),
                   std::move(emf)) {}

void BasicUpdator3D::updateE() {
  auto nx{_grid_space->sizeX()};
  auto ny{_grid_space->sizeY()};
  auto nz{_grid_space->sizeZ()};

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

  for (std::size_t i{0}; i < nx; ++i) {
    for (std::size_t j{1}; j < ny; ++j) {
      for (std::size_t k{1}; k < nz; ++k) {
        ex(i, j, k) = cexe(i, j, k) * ex(i, j, k) +
                      cexhy(i, j, k) * (hy(i, j, k) - hy(i, j, k - 1)) +
                      cexhz(i, j, k) * (hz(i, j, k) - hz(i, j - 1, k));
      }
    }
  }

  for (std::size_t i{1}; i < nx; ++i) {
    for (std::size_t j{0}; j < ny; ++j) {
      for (std::size_t k{1}; k < nz; ++k) {
        ey(i, j, k) = ceye(i, j, k) * ey(i, j, k) +
                      ceyhz(i, j, k) * (hz(i, j, k) - hz(i - 1, j, k)) +
                      ceyhx(i, j, k) * (hx(i, j, k) - hx(i, j, k - 1));
      }
    }
  }

  for (std::size_t i{1}; i < nx; ++i) {
    for (std::size_t j{1}; j < ny; ++j) {
      for (std::size_t k{0}; k < nz; ++k) {
        ez(i, j, k) = ceze(i, j, k) * ez(i, j, k) +
                      cezhx(i, j, k) * (hx(i, j, k) - hx(i, j - 1, k)) +
                      cezhy(i, j, k) * (hy(i, j, k) - hy(i - 1, j, k));
      }
    }
  }
}

}  // namespace xfdtd
