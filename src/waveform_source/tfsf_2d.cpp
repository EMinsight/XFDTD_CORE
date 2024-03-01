#include "xfdtd/waveform_source/tfsf_2d.h"

#include <memory>

#include "divider/divider.h"
#include "waveform_source/waveform_source_corrector.h"
#include "xfdtd/util/constant.h"

namespace xfdtd {

TFSF2D::TFSF2D(std::size_t distance_x, std::size_t distance_y, double phi,
               std::unique_ptr<Waveform> waveform)
    : TFSF{distance_x, distance_y,         0, constant::PI * 0.5, phi,
           0,          std::move(waveform)} {}

void TFSF2D::init(std::shared_ptr<GridSpace> grid_space,
                  std::shared_ptr<CalculationParam> calculation_param,
                  std::shared_ptr<EMF> emf) {
  defaultInit(std::move(grid_space), std::move(calculation_param),
              std::move(emf));
}

void TFSF2D::correctH() {
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
  for (auto j{js}; j < je + 1; ++j) {
    auto &hy_xn{emfPtr()->hy()(is - 1, j, 0)};
    const auto ez_i{ezInc(is, j, 0)};
    hy_xn -= cb_x * ez_i;
  }

  // xp
  for (auto j{js}; j < je + 1; ++j) {
    auto &hy_xp{emfPtr()->hy()(ie, j, 0)};
    const auto ez_i{ezInc(ie, j, 0)};
    hy_xp += cb_x * ez_i;
  }

  // yn
  for (auto i{is}; i < ie + 1; ++i) {
    auto &hx_yn{emfPtr()->hx()(i, js - 1, 0)};
    const auto ez_i{ezInc(i, js, 0)};
    hx_yn += cb_y * ez_i;
  }

  // yp
  for (auto i{is}; i < ie + 1; ++i) {
    auto &hx_yp{emfPtr()->hx()(i, je, 0)};
    const auto ez_i{ezInc(i, je, 0)};
    hx_yp -= cb_y * ez_i;
  }
}

void TFSF2D::correctE() {
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
  for (auto j{js}; j < je + 1; ++j) {
    auto &ez_xn{emfPtr()->ez()(is, j, 0)};
    const auto hy_i{hyInc(is - 1, j, 0)};
    ez_xn -= ca_x * hy_i;
  }

  // xp
  for (auto j{js}; j < je + 1; ++j) {
    auto &ez_xp{emfPtr()->ez()(ie, j, 0)};
    const auto hy_i{hyInc(ie, j, 0)};
    ez_xp += ca_x * hy_i;
  }

  // yn
  for (auto i{is}; i < ie + 1; ++i) {
    auto &ez_yn{emfPtr()->ez()(i, js, 0)};
    const auto hx_i{hxInc(i, js - 1, 0)};
    ez_yn += ca_y * hx_i;
  }

  // yp
  for (auto i{is}; i < ie + 1; ++i) {
    auto &ez_yp{emfPtr()->ez()(i, je, 0)};
    const auto hx_i{hxInc(i, je, 0)};
    ez_yp -= ca_y * hx_i;
  }
}

std::unique_ptr<Corrector> TFSF2D::generateCorrector(
    const Divider::IndexTask &task) {
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

  auto xn_intersection = Divider::taskIntersection(task, xn_task);
  auto xp_intersection = Divider::taskIntersection(task, xp_task);
  auto yn_intersection = Divider::taskIntersection(task, yn_task);
  auto yp_intersection = Divider::taskIntersection(task, yp_task);

  xn_task = Divider::makeTask(Divider::makeRange(is, is + 1),
                              Divider::makeRange<std::size_t>(1, 0),
                              Divider::makeRange<std::size_t>(1, 0));
  std::size_t j_s_xn;
  std::size_t j_e_xn;
  std::size_t k_s_xn;
  std::size_t k_e_xn;

  if (xn_intersection.has_value()) {
    j_s_xn = xn_intersection->_y_range[0];
    j_e_xn = xn_intersection->_y_range[1];
    k_s_xn = xn_intersection->_z_range[0];
    k_e_xn = xn_intersection->_z_range[1];
    xn_task = Divider::makeTask(Divider::makeRange(is, is + 1),
                                Divider::makeRange(j_s_xn, j_e_xn),
                                Divider::makeRange(k_s_xn, k_e_xn));
  }

  xp_task = Divider::makeTask(Divider::makeRange(ie, ie + 1),
                              Divider::makeRange<std::size_t>(1, 0),
                              Divider::makeRange<std::size_t>(1, 0));
  std::size_t j_s_xp;
  std::size_t j_e_xp;
  std::size_t k_s_xp;
  std::size_t k_e_xp;

  if (xp_intersection.has_value()) {
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

  auto zn_task = Divider::makeTask(Divider::makeRange<std::size_t>(1, 0),
                                   Divider::makeRange<std::size_t>(1, 0),
                                   Divider::makeRange(ks, ks + 1));

  auto zp_task = Divider::makeTask(Divider::makeRange<std::size_t>(1, 0),
                                   Divider::makeRange<std::size_t>(1, 0),
                                   Divider::makeRange(ks, ks + 1));

  auto domain_task = Divider::taskIntersection(
                         task, Divider::makeTask(Divider::makeRange(is, ie),
                                                 Divider::makeRange(js, je),
                                                 Divider::makeRange(ks, ke)))
                         .value();

  return std::make_unique<TFSF2DCorrector>(
      domain_task,
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
