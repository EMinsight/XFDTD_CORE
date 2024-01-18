#include "xfdtd/waveform_source/tfsf_2d.h"

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

}  // namespace xfdtd
