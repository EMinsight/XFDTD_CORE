#include "xfdtd/waveform_source/tfsf_3d.h"

#include <utility>

#include "waveform_source/waveform_source_corrector.h"

namespace xfdtd {

TFSF3D::TFSF3D(std::size_t x, std::size_t y, std::size_t z, double theta,
               double phi, double psi, std::unique_ptr<Waveform> waveform)
    : TFSF{x, y, z, theta, phi, psi, std::move(waveform)} {}

void TFSF3D::init(std::shared_ptr<GridSpace> grid_space,
                  std::shared_ptr<CalculationParam> calculation_param,
                  std::shared_ptr<EMF> emf) {
  defaultInit(std::move(grid_space), std::move(calculation_param),
              std::move(emf));
}

void TFSF3D::correctH() {
  const auto is{boxPtr()->origin().i()};
  const auto js{boxPtr()->origin().j()};
  const auto ks{boxPtr()->origin().k()};
  const auto ie{boxPtr()->end().i()};
  const auto je{boxPtr()->end().j()};
  const auto ke{boxPtr()->end().k()};

  auto emf{emfPtr()};
  const auto cb_x{cbx()};
  const auto cb_y{cby()};
  const auto cb_z{cbz()};

  // xn
  for (auto j{js}; j < je; ++j) {
    for (auto k{ks}; k < ke + 1; ++k) {
      auto& hz_xn{emfPtr()->hz()(is - 1, j, k)};
      const auto ey_i{eyInc(is, j, k)};
      hz_xn += cb_x * ey_i;
    }
  }
  for (auto j{js}; j < je + 1; ++j) {
    for (auto k{ks}; k < ke; ++k) {
      auto& hy_xn{emfPtr()->hy()(is - 1, j, k)};
      const auto ez_i{ezInc(is, j, k)};
      hy_xn -= cb_x * ez_i;
    }
  }
  // xp
  for (auto j{js}; j < je; ++j) {
    for (auto k{ks}; k < ke + 1; ++k) {
      auto& hz_xp{emfPtr()->hz()(ie, j, k)};
      const auto ey_i{eyInc(ie, j, k)};
      hz_xp -= cb_x * ey_i;
    }
  }
  for (auto j{js}; j < je + 1; ++j) {
    for (auto k{ks}; k < ke; ++k) {
      auto& hy_xp{emfPtr()->hy()(ie, j, k)};
      const auto ez_i{ezInc(ie, j, k)};
      hy_xp += cb_x * ez_i;
    }
  }

  // yn
  for (auto i{is}; i < ie + 1; ++i) {
    for (auto k{ks}; k < ke; ++k) {
      auto& hx_yn{emfPtr()->hx()(i, js - 1, k)};
      const auto ez_i{ezInc(i, js, k)};
      hx_yn += cb_y * ez_i;
    }
  }
  for (auto i{is}; i < ie; ++i) {
    for (auto k{ks}; k < ke + 1; ++k) {
      auto& hz_yn{emfPtr()->hz()(i, js - 1, k)};
      const auto ex_i{exInc(i, js, k)};
      hz_yn -= cb_y * ex_i;
    }
  }
  // yp
  for (auto i{is}; i < ie + 1; ++i) {
    for (auto k{ks}; k < ke; ++k) {
      auto& hx_yp{emfPtr()->hx()(i, je, k)};
      const auto ez_i{ezInc(i, je, k)};
      hx_yp -= cb_y * ez_i;
    }
  }
  for (auto i{is}; i < ie; ++i) {
    for (auto k{ks}; k < ke + 1; ++k) {
      auto& hz_yp{emfPtr()->hz()(i, je, k)};
      const auto ex_i{exInc(i, je, k)};
      hz_yp += cb_y * ex_i;
    }
  }

  // zn
  for (auto i{is}; i < ie; ++i) {
    for (auto j{js}; j < je + 1; ++j) {
      auto& hy_zn{emfPtr()->hy()(i, j, ks - 1)};
      const auto ex_i{exInc(i, j, ks)};
      hy_zn += cb_z * ex_i;
    }
  }
  for (auto i{is}; i < ie + 1; ++i) {
    for (auto j{js}; j < je; ++j) {
      auto& hx_zn{emfPtr()->hx()(i, j, ks - 1)};
      const auto ey_i{eyInc(i, j, ks)};
      hx_zn -= cb_z * ey_i;
    }
  }
  // zp
  for (auto i{is}; i < ie; ++i) {
    for (auto j{js}; j < je + 1; ++j) {
      auto& hy_zp{emfPtr()->hy()(i, j, ke)};
      const auto ex_i{exInc(i, j, ke)};
      hy_zp -= cb_z * ex_i;
    }
  }
  for (auto i{is}; i < ie + 1; ++i) {
    for (auto j{js}; j < je; ++j) {
      auto& hx_zp{emfPtr()->hx()(i, j, ke)};
      const auto ey_i{eyInc(i, j, ke)};
      hx_zp += cb_z * ey_i;
    }
  }
}

void TFSF3D::correctE() {
  const auto is{boxPtr()->origin().i()};
  const auto js{boxPtr()->origin().j()};
  const auto ks{boxPtr()->origin().k()};
  const auto ie{boxPtr()->end().i()};
  const auto je{boxPtr()->end().j()};
  const auto ke{boxPtr()->end().k()};

  auto emf{emfPtr()};
  const auto ca_x{cax()};
  const auto ca_y{cay()};
  const auto ca_z{caz()};

  // xn
  for (auto j{js}; j < je; ++j) {
    for (auto k{ks}; k < ke + 1; ++k) {
      auto& ey_xn{emfPtr()->ey()(is, j, k)};
      const auto hz_i{hzInc(is - 1, j, k)};
      ey_xn += ca_x * hz_i;
    }
  }
  for (auto j{js}; j < je + 1; ++j) {
    for (auto k{ks}; k < ke; ++k) {
      auto& ez_xn{emfPtr()->ez()(is, j, k)};
      const auto hy_i{hyInc(is - 1, j, k)};
      ez_xn -= ca_x * hy_i;
    }
  }
  // xp
  for (auto j{js}; j < je; ++j) {
    for (auto k{ks}; k < ke + 1; ++k) {
      auto& ey_xp{emfPtr()->ey()(ie, j, k)};
      const auto hz_i{hzInc(ie, j, k)};
      ey_xp -= ca_x * hz_i;
    }
  }
  for (auto j{js}; j < je + 1; ++j) {
    for (auto k{ks}; k < ke; ++k) {
      auto& ez_xp{emfPtr()->ez()(ie, j, k)};
      const auto hy_i{hyInc(ie, j, k)};
      ez_xp += ca_x * hy_i;
    }
  }

  // yn
  for (auto i{is}; i < ie + 1; ++i) {
    for (auto k{ks}; k < ke; ++k) {
      auto& ez_yn{emfPtr()->ez()(i, js, k)};
      const auto hx_i{hxInc(i, js - 1, k)};
      ez_yn += ca_y * hx_i;
    }
  }
  for (auto i{is}; i < ie; ++i) {
    for (auto k{ks}; k < ke + 1; ++k) {
      auto& ex_yn{emfPtr()->ex()(i, js, k)};
      const auto hz_i{hzInc(i, js - 1, k)};
      ex_yn -= ca_y * hz_i;
    }
  }
  // yp
  for (auto i{is}; i < ie + 1; ++i) {
    for (auto k{ks}; k < ke; ++k) {
      auto& ez_yp{emfPtr()->ez()(i, je, k)};
      const auto hx_i{hxInc(i, je, k)};
      ez_yp -= ca_y * hx_i;
    }
  }
  for (auto i{is}; i < ie; ++i) {
    for (auto k{ks}; k < ke + 1; ++k) {
      auto& ex_yp{emfPtr()->ex()(i, je, k)};
      const auto hz_i{hzInc(i, je, k)};
      ex_yp += ca_y * hz_i;
    }
  }

  // zn
  for (auto i{is}; i < ie; ++i) {
    for (auto j{js}; j < je + 1; ++j) {
      auto& ex_zn{emfPtr()->ex()(i, j, ks)};
      const auto hy_i{hyInc(i, j, ks - 1)};
      ex_zn += ca_z * hy_i;
    }
  }
  for (auto i{is}; i < ie + 1; ++i) {
    for (auto j{js}; j < je; ++j) {
      auto& ey_zn{emfPtr()->ey()(i, j, ks)};
      const auto hx_i{hxInc(i, j, ks - 1)};
      ey_zn -= ca_z * hx_i;
    }
  }
  // zp
  for (auto i{is}; i < ie; ++i) {
    for (auto j{js}; j < je + 1; ++j) {
      auto& ex_zp{emfPtr()->ex()(i, j, ke)};
      const auto hy_i{hyInc(i, j, ke)};
      ex_zp -= ca_z * hy_i;
    }
  }
  for (auto i{is}; i < ie + 1; ++i) {
    for (auto j{js}; j < je; ++j) {
      auto& ey_zp{emfPtr()->ey()(i, j, ke)};
      const auto hx_i{hxInc(i, j, ke)};
      ey_zp += ca_z * hx_i;
    }
  }
}

std::unique_ptr<Corrector> TFSF3D::generateCorrector(
    const Divider::IndexTask& task) {
  const auto is{boxPtr()->origin().i()};
  const auto js{boxPtr()->origin().j()};
  const auto ks{boxPtr()->origin().k()};
  const auto ie{boxPtr()->end().i()};
  const auto je{boxPtr()->end().j()};
  const auto ke{boxPtr()->end().k()};

  // xn
  auto xn_task =
      Divider::makeTask(Divider::makeRange(is, is + 1),
                        Divider::makeRange(js, je), Divider::makeRange(ks, ke));
  // xp
  auto xp_task =
      Divider::makeTask(Divider::makeRange(ie, ie + 1),
                        Divider::makeRange(js, je), Divider::makeRange(ks, ke));
  // yn
  auto yn_task = Divider::makeTask(Divider::makeRange(is, ie),
                                   Divider::makeRange(js, js + 1),
                                   Divider::makeRange(ks, ke));
  // yp
  auto yp_task = Divider::makeTask(Divider::makeRange(is, ie),
                                   Divider::makeRange(je, je + 1),
                                   Divider::makeRange(ks, ke));
  // zn
  auto zn_task =
      Divider::makeTask(Divider::makeRange(is, ie), Divider::makeRange(js, je),
                        Divider::makeRange(ks, ks + 1));
  // zp
  auto zp_task =
      Divider::makeTask(Divider::makeRange(is, ie), Divider::makeRange(js, je),
                        Divider::makeRange(ke, ke + 1));

  auto xn_intersection = Divider::taskIntersection(task, xn_task);
  auto xp_intersection = Divider::taskIntersection(task, xp_task);
  auto yn_intersection = Divider::taskIntersection(task, yn_task);
  auto yp_intersection = Divider::taskIntersection(task, yp_task);
  auto zn_intersection = Divider::taskIntersection(task, zn_task);
  auto zp_intersection = Divider::taskIntersection(task, zp_task);

  xn_task = Divider::makeTask(Divider::makeRange(is, is + 1),
                              Divider::makeRange<std::size_t>(1, 0),
                              Divider::makeRange<std::size_t>(1, 0));
  std::size_t j_s_xn;
  std::size_t j_e_xn;
  std::size_t k_s_xn;
  std::size_t k_e_xn;
  if (xn_intersection) {
    j_s_xn = xn_intersection->_y_range[0];
    j_e_xn = xn_intersection->_y_range[1];
    k_s_xn = xn_intersection->_z_range[0];
    k_e_xn = xn_intersection->_z_range[1];
    xn_task = Divider::makeTask(Divider::makeRange(is, is + 1),
                                Divider::makeRange(j_s_xn, j_e_xn),
                                Divider::makeRange(k_s_xn, k_e_xn));
  };

  xp_task = Divider::makeTask(Divider::makeRange(ie, ie + 1),
                              Divider::makeRange<std::size_t>(1, 0),
                              Divider::makeRange<std::size_t>(1, 0));

  std::size_t j_s_xp;
  std::size_t j_e_xp;
  std::size_t k_s_xp;
  std::size_t k_e_xp;
  if (xp_intersection) {
    j_s_xp = xp_intersection->_y_range[0];
    j_e_xp = xp_intersection->_y_range[1];
    k_s_xp = xp_intersection->_z_range[0];
    k_e_xp = xp_intersection->_z_range[1];
    xp_task = Divider::makeTask(Divider::makeRange(ie, ie + 1),
                                Divider::makeRange(j_s_xp, j_e_xp),
                                Divider::makeRange(k_s_xp, k_e_xp));
  }

  yn_task = Divider::makeTask(Divider::makeRange<std::size_t>(1, 0),
                              Divider::makeRange(js, js + 1),
                              Divider::makeRange<std::size_t>(1, 0));
  std::size_t i_s_yn;
  std::size_t i_e_yn;
  std::size_t k_s_yn;
  std::size_t k_e_yn;

  if (yn_intersection.has_value()) {
    i_s_yn = yn_intersection->_x_range[0];
    i_e_yn = yn_intersection->_x_range[1];
    k_s_yn = yn_intersection->_z_range[0];
    k_e_yn = yn_intersection->_z_range[1];
    yn_task = Divider::makeTask(Divider::makeRange(i_s_yn, i_e_yn),
                                Divider::makeRange(js, js + 1),
                                Divider::makeRange(k_s_yn, k_e_yn));
  }

  yp_task = Divider::makeTask(Divider::makeRange<std::size_t>(1, 0),
                              Divider::makeRange(je, je + 1),
                              Divider::makeRange<std::size_t>(1, 0));

  std::size_t i_s_yp;
  std::size_t i_e_yp;
  std::size_t k_s_yp;
  std::size_t k_e_yp;

  if (yp_intersection.has_value()) {
    i_s_yp = yp_intersection->_x_range[0];
    i_e_yp = yp_intersection->_x_range[1];
    k_s_yp = yp_intersection->_z_range[0];
    k_e_yp = yp_intersection->_z_range[1];
    yp_task = Divider::makeTask(Divider::makeRange(i_s_yp, i_e_yp),
                                Divider::makeRange(je, je + 1),
                                Divider::makeRange(k_s_yp, k_e_yp));
  }

  zn_task = Divider::makeTask(Divider::makeRange<std::size_t>(1, 0),
                              Divider::makeRange<std::size_t>(1, 0),
                              Divider::makeRange(ks, ks + 1));
  std::size_t i_s_zn;
  std::size_t i_e_zn;
  std::size_t j_s_zn;
  std::size_t j_e_zn;
  if (zn_intersection.has_value()) {
    i_s_zn = zn_intersection->_x_range[0];
    i_e_zn = zn_intersection->_x_range[1];
    j_s_zn = zn_intersection->_y_range[0];
    j_e_zn = zn_intersection->_y_range[1];
    zn_task = Divider::makeTask(Divider::makeRange(i_s_zn, i_e_zn),
                                Divider::makeRange(j_s_zn, j_e_zn),
                                Divider::makeRange(ks, ks + 1));
  }

  zp_task = Divider::makeTask(Divider::makeRange<std::size_t>(1, 0),
                              Divider::makeRange<std::size_t>(1, 0),
                              Divider::makeRange(ke, ke + 1));
  std::size_t i_s_zp;
  std::size_t i_e_zp;
  std::size_t j_s_zp;
  std::size_t j_e_zp;
  if (zp_intersection.has_value()) {
    i_s_zp = zp_intersection->_x_range[0];
    i_e_zp = zp_intersection->_x_range[1];
    j_s_zp = zp_intersection->_y_range[0];
    j_e_zp = zp_intersection->_y_range[1];
    zp_task = Divider::makeTask(Divider::makeRange(i_s_zp, i_e_zp),
                                Divider::makeRange(j_s_zp, j_e_zp),
                                Divider::makeRange(ke, ke + 1));
  }

  auto domain_task = Divider::taskIntersection(
      task,
      Divider::makeTask(Divider::makeRange(is, ie), Divider::makeRange(js, je),
                        Divider::makeRange(ks, ke)));

  return std::make_unique<TFSF3DCorrector>(
      domain_task.value(),
      Divider::makeTask(Divider::makeRange<std::size_t>(0, 0),
                        Divider::makeRange<std::size_t>(0, 0),
                        Divider::makeRange<std::size_t>(0, 0)),
      gridSpace(), calculationParam(), emf(), waveform()->value(),
      boxPtr()->origin(), std::move(xn_task), std::move(xp_task),
      std::move(yn_task), std::move(yp_task), std::move(zn_task),
      std::move(zp_task), _projection_x_int, _projection_y_int,
      _projection_z_int, _projection_x_half, _projection_y_half,
      _projection_z_half, _ex_inc, _ey_inc, _ez_inc, _hx_inc, _hy_inc, _hz_inc,
      cax(), cbx(), cay(), cby(), caz(), cbz());
}

}  // namespace xfdtd
