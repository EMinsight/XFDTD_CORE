#include "updator/basic_updator.h"
#include "updator/update_scheme.h"

namespace xfdtd {

BasicUpdator3D::BasicUpdator3D(
    std::shared_ptr<const GridSpace> grid_space,
    std::shared_ptr<const CalculationParam> calculation_param,
    std::shared_ptr<EMF> emf, Divider::IndexTask task)
    : BasicUpdator(std::move(grid_space), std::move(calculation_param),
                   std::move(emf), task) {}

std::string BasicUpdator3D::toString() const {
  std::stringstream ss;
  ss << BasicUpdator::toString() << "\n";
  ss << "E field: \n";
  auto x_range = task()._x_range;
  auto y_range = task()._y_range;
  auto z_range = task()._z_range;

  auto ex_task =
      Divider::makeTask(Divider::makeRange(x_range.start(), x_range.end()),
                        Divider::makeRange(y_range.start() + 1, y_range.end()),
                        Divider::makeRange(z_range.start() + 1, z_range.end()));
  ss << "Ex: " << ex_task.toString() << "\n";
  auto ey_task =
      Divider::makeTask(Divider::makeRange(x_range.start() + 1, x_range.end()),
                        Divider::makeRange(y_range.start(), y_range.end()),
                        Divider::makeRange(z_range.start() + 1, z_range.end()));
  ss << "Ey: " << ey_task.toString() << "\n";
  auto ez_task =
      Divider::makeTask(Divider::makeRange(x_range.start() + 1, x_range.end()),
                        Divider::makeRange(y_range.start() + 1, y_range.end()),
                        Divider::makeRange(z_range.start(), z_range.end()));
  ss << "Ez: " << ez_task.toString();
  return ss.str();
}

void BasicUpdator3D::updateE() {
  const auto is = task()._x_range[0];
  const auto ie = task()._x_range[1];
  const auto js = task()._y_range[0];
  const auto je = task()._y_range[1];
  const auto ks = task()._z_range[0];
  const auto ke = task()._z_range[1];

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

  for (std::size_t i{is}; i < ie; ++i) {
    for (std::size_t j{js + 1}; j < je; ++j) {
      for (std::size_t k{ks + 1}; k < ke; ++k) {
        ex(i, j, k) = eNext(cexe(i, j, k), ex(i, j, k), cexhy(i, j, k),
                            hy(i, j, k), hy(i, j, k - 1), cexhz(i, j, k),
                            hz(i, j, k), hz(i, j - 1, k));
      }
    }
  }

  for (std::size_t i{is + 1}; i < ie; ++i) {
    for (std::size_t j{js}; j < je; ++j) {
      for (std::size_t k{ks + 1}; k < ke; ++k) {
        ey(i, j, k) = eNext(ceye(i, j, k), ey(i, j, k), ceyhz(i, j, k),
                            hz(i, j, k), hz(i - 1, j, k), ceyhx(i, j, k),
                            hx(i, j, k), hx(i, j, k - 1));
      }
    }
  }

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
    const auto& hy_buffer = hyBufferZN();
    const auto& hz_buffer = hzBufferYN();
    updateExLineX(is, ie, js, ks, cexe, cexhy, cexhz, hy, hy_buffer, hz,
                  hz_buffer, ex);
  }

  if (!contain_xn_edge && !contain_zn_edge) {
    const auto& hz_buffer = hzBufferXN();
    const auto& hx_buffer = hxBufferZN();
    updateEyLineY(is, js, je, ks, ceye, ceyhz, ceyhx, hz, hz_buffer, hx,
                  hx_buffer, ey);
  }

  if (!contain_xn_edge && !contain_yn_edge) {
    const auto& hx_buffer = hxBufferYN();
    const auto& hy_buffer = hyBufferXN();
    updateEzLineZ(is, js, ke, ks, ceze, cezhx, cezhy, hx, hx_buffer, hy,
                  hy_buffer, ez);
  }

  if (!contain_xn_edge) {
    const auto& hy_buffer = hyBufferXN();
    const auto& hz_buffer = hzBufferXN();
    updateEyEdgeYZ(is, js, je, ks, ke, ceye, ceyhz, ceyhx, hz, hz_buffer, hx,
                   ey);
    updateEzEdgeYZ(is, js, je, ks, ke, ceze, cezhx, cezhy, hx, hy, hy_buffer,
                   ez);
  }

  if (!contain_yn_edge) {
    const auto& hz_buffer = hzBufferYN();
    const auto& hx_buffer = hxBufferYN();
    updateEzEdgeXZ(is, ie, js, ks, ke, ceze, cezhx, cezhy, hx, hx_buffer, hy,
                   ez);
    updateExEdgeXZ(is, ie, js, ks, ke, cexe, cexhy, cexhz, hy, hz, hz_buffer,
                   ex);
  }

  if (!contain_zn_edge) {
    const auto& hy_buffer = hyBufferZN();
    const auto& hx_buffer = hxBufferZN();
    updateExEdgeXY(is, ie, js, je, ks, cexe, cexhy, cexhz, hy, hy_buffer, hz,
                   ex);
    updateEyEdgeXY(is, ie, js, je, ks, ceye, ceyhz, ceyhx, hz, hx, hx_buffer,
                   ey);
  }
}

}  // namespace xfdtd
