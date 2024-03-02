#include <sstream>
#include <utility>
#include <xtensor.hpp>

#include "divider/divider.h"
#include "updator/basic_updator.h"
#include "updator/update_scheme.h"
#include "updator/updator.h"

namespace xfdtd {

BasicUpdatorTE::BasicUpdatorTE(
    std::shared_ptr<const GridSpace> grid_space,
    std::shared_ptr<const CalculationParam> calculation_param,
    std::shared_ptr<EMF> emf, Divider::IndexTask task)
    : BasicUpdator(std::move(grid_space), std::move(calculation_param),
                   std::move(emf), task) {}

void BasicUpdatorTE::updateE() {
  const auto is = task()._x_range[0];
  const auto ie = task()._x_range[1];
  const auto js = task()._y_range[0];
  const auto je = task()._y_range[1];
  const auto ks = task()._z_range[0];
  const auto ke = task()._z_range[1];

  const auto& ceze{_calculation_param->fdtdCoefficient()->ceze()};
  const auto& cezhx{_calculation_param->fdtdCoefficient()->cezhx()};
  const auto& cezhy{_calculation_param->fdtdCoefficient()->cezhy()};

  auto& hx{_emf->hx()};
  auto& hy{_emf->hy()};
  auto& ez{_emf->ez()};

  for (std::size_t i{is + 1}; i < ie; ++i) {
    for (std::size_t j{js + 1}; j < je; ++j) {
      for (std::size_t k{ks}; k < ke; ++k) {
        ez(i, j, k) = eNext(ceze(i, j, k), ez(i, j, k), cezhx(i, j, k),
                            hx(i, j, k), hx(i, j - 1, k), cezhy(i, j, k),
                            hy(i, j, k), hy(i - 1, j, k));
      }
    }
  }

  updateEEdge();
}

std::string BasicUpdatorTE::toString() const {
  std::stringstream ss;
  ss << Updator::toString() << "\n";
  auto task_ez = Divider::makeTask<std::size_t>(
      {task().xRange().start() + 1, task().xRange().end()},
      {task().yRange().start() + 1, task().yRange().end()},
      {task().zRange().start(), task().zRange().end()});
  ss << "Ez field: " << task_ez.toString();
  return ss.str();
}

void BasicUpdatorTE::updateEEdge() {
  const auto is = task()._x_range[0];
  const auto ie = task()._x_range[1];
  const auto js = task()._y_range[0];
  const auto je = task()._y_range[1];
  const auto ks = task()._z_range[0];
  const auto ke = task()._z_range[1];

  const auto& ceze{_calculation_param->fdtdCoefficient()->ceze()};
  const auto& cezhx{_calculation_param->fdtdCoefficient()->cezhx()};
  const auto& cezhy{_calculation_param->fdtdCoefficient()->cezhy()};

  const auto& hx{_emf->hx()};
  const auto& hy{_emf->hy()};
  auto& ez{_emf->ez()};

  const auto& hy_buffer = hyBuffer();
  const auto& hx_buffer = hxBuffer();

  bool contain_xn_edge = containXNEdge();
  bool contain_yn_edge = containYNEdge();

  if (!contain_xn_edge && !contain_yn_edge) {
    updateEzCornerXY(is, js, ks, ke, ceze, cezhx, cezhy, hx, hx_buffer, hy,
                     hy_buffer, ez);
  }

  if (!contain_xn_edge) {
    updateEzEdgeYZ(is, js, je, ks, ke, ceze, cezhx, cezhy, hx, hy, hy_buffer,
                   ez);
  }

  if (!contain_yn_edge) {
    updateEzEdgeXZ(is, ie, js, ks, ke, ceze, cezhx, cezhy, hx, hx_buffer, hy,
                   ez);
  }
}

}  // namespace xfdtd
