#include <xfdtd/util/fdtd_basic.h>

#include "updator/basic_updator.h"
#include "updator/update_scheme.h"

namespace xfdtd {

BasicUpdator3D::BasicUpdator3D(
    std::shared_ptr<const GridSpace> grid_space,
    std::shared_ptr<const CalculationParam> calculation_param,
    std::shared_ptr<EMF> emf, IndexTask task)
    : BasicUpdator(std::move(grid_space), std::move(calculation_param),
                   std::move(emf), task) {}

std::string BasicUpdator3D::toString() const {
  std::stringstream ss;
  ss << BasicUpdator::toString() << "\n";
  ss << "E field: \n";
  auto x_range = task()._x_range;
  auto y_range = task()._y_range;
  auto z_range = task()._z_range;

  auto ex_task = makeTask(makeRange(x_range.start(), x_range.end()),
                          makeRange(y_range.start() + 1, y_range.end()),
                          makeRange(z_range.start() + 1, z_range.end()));
  ss << "Ex: " << ex_task.toString() << "\n";
  auto ey_task = makeTask(makeRange(x_range.start() + 1, x_range.end()),
                          makeRange(y_range.start(), y_range.end()),
                          makeRange(z_range.start() + 1, z_range.end()));
  ss << "Ey: " << ey_task.toString() << "\n";
  auto ez_task = makeTask(makeRange(x_range.start() + 1, x_range.end()),
                          makeRange(y_range.start() + 1, y_range.end()),
                          makeRange(z_range.start(), z_range.end()));
  ss << "Ez: " << ez_task.toString();
  return ss.str();
}

void BasicUpdator3D::updateE() {
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
  // auto js = basic::GridStructure::exFDTDUpdateYStart(task.yRange().start());
  auto js = task.yRange().start() == 0 ? 1 : task.yRange().start();
  auto je = basic::GridStructure::exFDTDUpdateYEnd(task.yRange().end());
  // auto ks = basic::GridStructure::exFDTDUpdateZStart(task.zRange().start());
  auto ks = task.zRange().start() == 0 ? 1 : task.zRange().start();
  auto ke = basic::GridStructure::exFDTDUpdateZEnd(task.zRange().end());
  for (std::size_t i{is}; i < ie; ++i) {
    for (std::size_t j{js}; j < je; ++j) {
      for (std::size_t k{ks}; k < ke; ++k) {
        ex(i, j, k) = eNext(cexe(i, j, k), ex(i, j, k), cexhy(i, j, k),
                            hy(i, j, k), hy(i, j, k - 1), cexhz(i, j, k),
                            hz(i, j, k), hz(i, j - 1, k));
      }
    }
  }

  // is = basic::GridStructure::eyFDTDUpdateXStart(task.xRange().start());
  is = task.xRange().start() == 0 ? 1 : task.xRange().start();
  ie = basic::GridStructure::eyFDTDUpdateXEnd(task.xRange().end());
  js = basic::GridStructure::eyFDTDUpdateYStart(task.yRange().start());
  je = basic::GridStructure::eyFDTDUpdateYEnd(task.yRange().end());
  // ks = basic::GridStructure::eyFDTDUpdateZStart(task.zRange().start());
  ks = task.zRange().start() == 0 ? 1 : task.zRange().start();
  ke = basic::GridStructure::eyFDTDUpdateZEnd(task.zRange().end());
  for (std::size_t i{is}; i < ie; ++i) {
    for (std::size_t j{js}; j < je; ++j) {
      for (std::size_t k{ks}; k < ke; ++k) {
        ey(i, j, k) = eNext(ceye(i, j, k), ey(i, j, k), ceyhz(i, j, k),
                            hz(i, j, k), hz(i - 1, j, k), ceyhx(i, j, k),
                            hx(i, j, k), hx(i, j, k - 1));
      }
    }
  }

  // is = basic::GridStructure::ezFDTDUpdateXStart(task.xRange().start());
  is = task.xRange().start() == 0 ? 1 : task.xRange().start();
  ie = basic::GridStructure::ezFDTDUpdateXEnd(task.xRange().end());
  // js = basic::GridStructure::ezFDTDUpdateYStart(task.yRange().start());
  js = task.yRange().start() == 0 ? 1 : task.yRange().start();
  je = basic::GridStructure::ezFDTDUpdateYEnd(task.yRange().end());
  ks = basic::GridStructure::ezFDTDUpdateZStart(task.zRange().start());
  ke = basic::GridStructure::ezFDTDUpdateZEnd(task.zRange().end());
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

void BasicUpdator3D::updateEEdge() {
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
      ex(i, j, k) =
          eNext(cexe(i, j, k), ex(i, j, k), cexhy(i, j, k), hy(i, j, k),
                hy(i, j, k - 1), cexhz(i, j, k), hz(i, j, k), hz(i, j - 1, k));
    }
  }

  if (!contain_xn_edge && !contain_zn_edge) {
    auto i = is;
    auto k = ks;
    for (std::size_t j{js}; j < je; ++j) {
      ey(i, j, k) =
          eNext(ceye(i, j, k), ey(i, j, k), ceyhz(i, j, k), hz(i, j, k),
                hz(i - 1, j, k), ceyhx(i, j, k), hx(i, j, k), hx(i, j, k - 1));
    }
  }

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
    for (std::size_t j{js}; j < je; ++j) {
      for (std::size_t k{ks + 1}; k < ke; ++k) {
        ey(i, j, k) = eNext(ceye(i, j, k), ey(i, j, k), ceyhz(i, j, k),
                            hz(i, j, k), hz(i - 1, j, k), ceyhx(i, j, k),
                            hx(i, j, k), hx(i, j, k - 1));
      }
    }
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
    for (std::size_t i{is}; i < ie; ++i) {
      for (std::size_t k{ks + 1}; k < ke; ++k) {
        ex(i, j, k) = eNext(cexe(i, j, k), ex(i, j, k), cexhy(i, j, k),
                            hy(i, j, k), hy(i, j, k - 1), cexhz(i, j, k),
                            hz(i, j, k), hz(i, j - 1, k));
      }
    }
  }

  if (!contain_zn_edge) {
    auto k = ks;
    for (std::size_t i{is}; i < ie; ++i) {
      for (std::size_t j{js + 1}; j < je; ++j) {
        ex(i, j, k) = eNext(cexe(i, j, k), ex(i, j, k), cexhy(i, j, k),
                            hy(i, j, k), hy(i, j, k - 1), cexhz(i, j, k),
                            hz(i, j, k), hz(i, j - 1, k));
      }
    }
    for (std::size_t i{is + 1}; i < ie; ++i) {
      for (std::size_t j{js}; j < je; ++j) {
        ey(i, j, k) = eNext(ceye(i, j, k), ey(i, j, k), ceyhz(i, j, k),
                            hz(i, j, k), hz(i - 1, j, k), ceyhx(i, j, k),
                            hx(i, j, k), hx(i, j, k - 1));
      }
    }
  }
}

}  // namespace xfdtd
