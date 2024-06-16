
#include <xfdtd/util/fdtd_basic.h>

#include <sstream>
#include <utility>
#include <xtensor.hpp>

#include "updator/basic_updator.h"
#include "updator/update_scheme.h"
#include "updator/updator.h"

namespace xfdtd {

BasicUpdatorTE::BasicUpdatorTE(
    std::shared_ptr<const GridSpace> grid_space,
    std::shared_ptr<const CalculationParam> calculation_param,
    std::shared_ptr<EMF> emf, IndexTask task)
    : BasicUpdator(std::move(grid_space), std::move(calculation_param),
                   std::move(emf), task) {}

void BasicUpdatorTE::updateE() {
  const auto task = this->task();
  const auto x_range = task.xRange();
  const auto y_range = task.yRange();
  const auto z_range = task.zRange();
  // const auto is = basic::GridStructure::ezFDTDUpdateXStart(x_range.start());
  const auto is = x_range.start() == 0 ? 1 : x_range.start();
  const auto ie = basic::GridStructure::ezFDTDUpdateXEnd(x_range.end());
  // const auto js = basic::GridStructure::ezFDTDUpdateYStart(y_range.start());
  const auto js = y_range.start() == 0 ? 1 : y_range.start();
  const auto je = basic::GridStructure::ezFDTDUpdateYEnd(y_range.end());
  const auto ks = basic::GridStructure::ezFDTDUpdateZStart(z_range.start());
  const auto ke = basic::GridStructure::ezFDTDUpdateZEnd(z_range.end());

  const auto& ceze{_calculation_param->fdtdCoefficient()->ceze()};
  const auto& cezhx{_calculation_param->fdtdCoefficient()->cezhx()};
  const auto& cezhy{_calculation_param->fdtdCoefficient()->cezhy()};

  auto& hx{_emf->hx()};
  auto& hy{_emf->hy()};
  auto& ez{_emf->ez()};

  for (std::size_t i{is}; i < ie; ++i) {
    for (std::size_t j{js}; j < je; ++j) {
      for (std::size_t k{ks}; k < ke; ++k) {
        ez(i, j, k) = eNext(ceze(i, j, k), ez(i, j, k), cezhx(i, j, k),
                            hx(i, j, k), hx(i, j - 1, k), cezhy(i, j, k),
                            hy(i, j, k), hy(i - 1, j, k));
      }
    }
  }

  // updateEEdge();
}

std::string BasicUpdatorTE::toString() const {
  std::stringstream ss;
  ss << Updator::toString() << "\n";
  auto task_ez = basic::GridStructure::ezFDTDUpdateTask(task());
  ss << " Ez field: " << task_ez.toString();
  return ss.str();
}

void BasicUpdatorTE::updateEEdge() {
  const auto is = task().xRange().start();
  const auto ie = task().xRange().end();
  const auto js = task().yRange().start();
  const auto je = task().yRange().end();
  const auto ks = task().zRange().start();
  const auto ke = task().zRange().end();

  const auto& ceze{_calculation_param->fdtdCoefficient()->ceze()};
  const auto& cezhx{_calculation_param->fdtdCoefficient()->cezhx()};
  const auto& cezhy{_calculation_param->fdtdCoefficient()->cezhy()};

  const auto& hx{_emf->hx()};
  const auto& hy{_emf->hy()};
  auto& ez{_emf->ez()};

  bool contain_xn_edge = containXNEdge();
  bool contain_yn_edge = containYNEdge();

  if (!contain_xn_edge && !contain_yn_edge) {
    auto i = is;
    auto j = js;
    for (std::size_t k{ks}; k < ke; ++k) {
      ez(i, j, k) =
          eNext(ceze(i, j, k), ez(i, j, k), cezhx(i, j, k), hx(i, j, k),
                hx(i, j - 1, k), cezhy(i, j, k), hy(i, j, k), hy(i - 1, j, k));
    }
  }

  if (!contain_xn_edge) {
    auto i = is;
    for (std::size_t j{js + 1}; j < je; ++j) {
      for (std::size_t k{ks}; k < ke; ++k) {
        ez(i, j, k) = eNext(ceze(i, j, k), ez(i, j, k), cezhx(i, j, k),
                            hx(i, j, k), hx(i, j - 1, k), cezhy(i, j, k),
                            hy(i, j, k), hy(i - 1, j, k));
      }
    }
  }

  if (!contain_yn_edge) {
    auto j = js;
    for (std::size_t i{is + 1}; i < ie; ++i) {
      for (std::size_t k{ks}; k < ke; ++k) {
        ez(i, j, k) = eNext(ceze(i, j, k), ez(i, j, k), cezhx(i, j, k),
                            hx(i, j, k), hx(i, j - 1, k), cezhy(i, j, k),
                            hy(i, j, k), hy(i - 1, j, k));
      }
    }
  }
}

}  // namespace xfdtd
