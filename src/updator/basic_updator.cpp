#include "updator/basic_updator.h"

#include <utility>
#include <xtensor.hpp>

#include "updator/update_scheme.h"
#include <xfdtd/divider/divider.h>
#include "updator/updator.h"

namespace xfdtd {

BasicUpdator::BasicUpdator(
    std::shared_ptr<const GridSpace> grid_space,
    std::shared_ptr<const CalculationParam> calculation_param,
    std::shared_ptr<EMF> emf, Divider::IndexTask task)
    : Updator(std::move(grid_space), std::move(calculation_param),
              std::move(emf), task) {}

void BasicUpdator::updateH() {
  const auto is = task()._x_range[0];
  const auto ie = task()._x_range[1];
  const auto js = task()._y_range[0];
  const auto je = task()._y_range[1];
  const auto ks = task()._z_range[0];
  const auto ke = task()._z_range[1];

  const auto& chxh{_calculation_param->fdtdCoefficient()->chxh()};
  const auto& chyh{_calculation_param->fdtdCoefficient()->chyh()};
  const auto& chzh{_calculation_param->fdtdCoefficient()->chzh()};
  const auto& chxey{_calculation_param->fdtdCoefficient()->chxey()};
  const auto& chxez{_calculation_param->fdtdCoefficient()->chxez()};
  const auto& chyex{_calculation_param->fdtdCoefficient()->chyex()};
  const auto& chyez{_calculation_param->fdtdCoefficient()->chyez()};
  const auto& chzex{_calculation_param->fdtdCoefficient()->chzex()};
  const auto& chzey{_calculation_param->fdtdCoefficient()->chzey()};

  auto& hx{_emf->hx()};
  auto& hy{_emf->hy()};
  auto& hz{_emf->hz()};
  const auto& ex{_emf->ex()};
  const auto& ey{_emf->ey()};
  const auto& ez{_emf->ez()};

  for (std::size_t i{is}; i < ie; ++i) {
    for (std::size_t j{js}; j < je; ++j) {
      for (std::size_t k{ks}; k < ke; ++k) {
        hx(i, j, k) =
            hNext(chxh(i, j, k), hx(i, j, k), chxey(i, j, k), ey(i, j, k + 1),
                  ey(i, j, k), chxez(i, j, k), ez(i, j + 1, k), ez(i, j, k));

        hy(i, j, k) =
            hNext(chyh(i, j, k), hy(i, j, k), chyez(i, j, k), ez(i + 1, j, k),
                  ez(i, j, k), chyex(i, j, k), ex(i, j, k + 1), ex(i, j, k));

        hz(i, j, k) =
            hNext(chzh(i, j, k), hz(i, j, k), chzex(i, j, k), ex(i, j + 1, k),
                  ex(i, j, k), chzey(i, j, k), ey(i + 1, j, k), ey(i, j, k));
      }
    }
  }

  updateHEdge();
}

BasicUpdatorTEM::BasicUpdatorTEM(
    std::shared_ptr<const GridSpace> grid_space,
    std::shared_ptr<const CalculationParam> calculation_param,
    std::shared_ptr<EMF> emf, Divider::IndexTask task)
    : BasicUpdator(std::move(grid_space), std::move(calculation_param),
                   std::move(emf), task) {}

void BasicUpdatorTEM::updateE() {
  const auto ks = task()._z_range[0];
  const auto ke = task()._z_range[1];

  const auto& cexe{_calculation_param->fdtdCoefficient()->cexe()};
  const auto& cexhy{_calculation_param->fdtdCoefficient()->cexhy()};

  auto& hy{_emf->hy()};
  auto& ex{_emf->ex()};

  for (std::size_t k{ks}; k < ke; ++k) {
    ex(0, 0, k) = eNext(cexe(0, 0, k), ex(0, 0, k), cexhy(0, 0, k), hy(0, 0, k),
                        hy(0, 0, k - 1), 0.0, 0.0, 0.0);
  }

  updateEEdge();
}

}  // namespace xfdtd
