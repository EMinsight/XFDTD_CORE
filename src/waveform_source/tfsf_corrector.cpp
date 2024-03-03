#include <memory>

#include "divider/divider.h"
#include "waveform_source/waveform_source_corrector.h"

namespace xfdtd {

double TFSFCorrector::exInc(std::size_t i, std::size_t j, std::size_t k) const {
  i = i - start().i() + 1;
  j = j - start().j();
  k = k - start().k();
  auto projection{_projection_x_half(i) + _projection_y_int(j) +
                  _projection_z_int(k)};
  auto index{static_cast<std::size_t>(projection)};
  auto weight{projection - index};
  return (1 - weight) * _ex_inc(index) + weight * _ex_inc(index + 1);
}

double TFSFCorrector::eyInc(std::size_t i, std::size_t j, std::size_t k) const {
  i = i - start().i();
  j = j - start().j() + 1;
  k = k - start().k();
  auto projection{_projection_x_int(i) + _projection_y_half(j) +
                  _projection_z_int(k)};
  auto index{static_cast<std::size_t>(projection)};
  auto weight{projection - index};
  return (1 - weight) * _ey_inc(index) + weight * _ey_inc(index + 1);
}

double TFSFCorrector::ezInc(std::size_t i, std::size_t j, std::size_t k) const {
  i = i - start().i();
  j = j - start().j();
  k = k - start().k() + 1;
  auto projection{_projection_x_int(i) + _projection_y_int(j) +
                  _projection_z_half(k)};
  auto index{static_cast<std::size_t>(projection)};
  auto weight{projection - index};
  return (1 - weight) * _ez_inc(index) + weight * _ez_inc(index + 1);
}

double TFSFCorrector::hxInc(std::size_t i, std::size_t j, std::size_t k) const {
  i = i - start().i();
  j = j - start().j() + 1;
  k = k - start().k() + 1;
  auto projection{_projection_x_int(i) + _projection_y_half(j) +
                  _projection_z_half(k) - 0.5};
  auto index{static_cast<std::size_t>(projection)};
  auto weight{projection - index};
  return (1 - weight) * _hx_inc(index) + weight * _hx_inc(index + 1);
}

double TFSFCorrector::hyInc(std::size_t i, std::size_t j, std::size_t k) const {
  i = i - start().i() + 1;
  j = j - start().j();
  k = k - start().k() + 1;
  auto projection{_projection_x_half(i) + _projection_y_int(j) +
                  _projection_z_half(k) - 0.5};
  auto index{static_cast<std::size_t>(projection)};
  auto weight{projection - index};
  return (1 - weight) * _hy_inc(index) + weight * _hy_inc(index + 1);
}

double TFSFCorrector::hzInc(std::size_t i, std::size_t j, std::size_t k) const {
  i = i - start().i() + 1;
  j = j - start().j() + 1;
  k = k - start().k();
  auto projection{_projection_x_half(i) + _projection_y_half(j) +
                  _projection_z_int(k) - 0.5};
  auto index{static_cast<std::size_t>(projection)};
  auto weight{projection - index};
  return (1 - weight) * _hz_inc(index) + weight * _hz_inc(index + 1);
}

void TFSF1DCorrector::correctE() {
  const auto lk = task()._z_range[0];
  const auto rk = task()._z_range[1];

  auto& ex_zn{emfPtr()->ex()(0, 0, lk)};
  auto hy_i{hyInc(0, 0, lk - 1)};
  ex_zn += cax() * hy_i;

  auto& ex_zp{emfPtr()->ex()(0, 0, rk)};
  hy_i = {hyInc(0, 0, rk)};
  ex_zp -= cax() * hy_i;
}

void TFSF1DCorrector::correctH() {
  const auto lk = task()._z_range[0];
  const auto rk = task()._z_range[1];

  auto& hy_zn{emfPtr()->hy()(0, 0, lk - 1)};
  auto ex_i{exInc(0, 0, lk)};
  hy_zn += cbx() * ex_i;

  auto& hy_zp{emfPtr()->hy()(0, 0, rk)};
  ex_i = exInc(0, 0, rk);
  hy_zp -= cbx() * ex_i;
}

void TFSF2DCorrector::correctE() {
  const auto is{task()._x_range[0]};
  const auto js{task()._y_range[0]};
  const auto ie{task()._x_range[1]};
  const auto je{task()._y_range[1]};

  auto emf{emfPtr()};
  const auto ca_x{cax()};
  const auto ca_y{cay()};

  // xn
  const auto js_xn = jStartXN();
  const auto je_xn = jEndXN();
  for (auto j{js_xn}; j < je_xn + 1; ++j) {
    auto& ez_xn{emfPtr()->ez()(is, j, 0)};
    const auto hy_i{hyInc(is - 1, j, 0)};
    ez_xn -= ca_x * hy_i;
  }

  // xp
  const auto js_xp = jStartXP();
  const auto je_xp = jEndXP();
  for (auto j{js_xp}; j < je_xp + 1; ++j) {
    auto& ez_xp{emfPtr()->ez()(ie, j, 0)};
    const auto hy_i{hyInc(ie, j, 0)};
    ez_xp += ca_x * hy_i;
  }

  // yn
  const auto is_yn = iStartYN();
  const auto ie_yn = iEndYN();
  for (auto i{is_yn}; i < ie_yn + 1; ++i) {
    auto& ez_yn{emfPtr()->ez()(i, js, 0)};
    const auto hx_i{hxInc(i, js - 1, 0)};
    ez_yn += ca_y * hx_i;
  }

  // yp
  const auto is_yp = iStartYP();
  const auto ie_yp = iEndYP();
  for (auto i{is_yp}; i < ie_yp + 1; ++i) {
    auto& ez_yp{emfPtr()->ez()(i, je, 0)};
    const auto hx_i{hxInc(i, je, 0)};
    ez_yp -= ca_y * hx_i;
  }
}

void TFSF2DCorrector::correctH() {
  const auto is{task()._x_range[0]};
  const auto js{task()._y_range[0]};
  const auto ie{task()._x_range[1]};
  const auto je{task()._y_range[1]};

  auto emf{emfPtr()};
  const auto cb_x{cbx()};
  const auto cb_y{cby()};

  // xn
  const auto js_xn = jStartXN();
  const auto je_xn = jEndXN();
  for (auto j{js_xn}; j < je_xn + 1; ++j) {
    auto& hy_xn{emfPtr()->hy()(is - 1, j, 0)};
    const auto ez_i{ezInc(is, j, 0)};
    hy_xn -= cb_x * ez_i;
  }

  // xp
  const auto js_xp = jStartXP();
  const auto je_xp = jEndXP();
  for (auto j{js_xp}; j < je_xp + 1; ++j) {
    auto& hy_xp{emfPtr()->hy()(ie, j, 0)};
    const auto ez_i{ezInc(ie, j, 0)};
    hy_xp += cb_x * ez_i;
  }

  // yn
  const auto is_yn = iStartYN();
  const auto ie_yn = iEndYN();
  for (auto i{is_yn}; i < ie_yn + 1; ++i) {
    auto& hx_yn{emfPtr()->hx()(i, js - 1, 0)};
    const auto ez_i{ezInc(i, js, 0)};
    hx_yn += cb_y * ez_i;
  }

  // yp
  const auto is_yp = iStartYP();
  const auto ie_yp = iEndYP();
  for (auto i{is_yp}; i < ie_yp + 1; ++i) {
    auto& hx_yp{emfPtr()->hx()(i, je, 0)};
    const auto ez_i{ezInc(i, je, 0)};
    hx_yp -= cb_y * ez_i;
  }
}

void TFSF3DCorrector::correctE() {
  const auto is{task()._x_range[0]};
  const auto js{task()._y_range[0]};
  const auto ks{task()._z_range[0]};
  const auto ie{task()._x_range[1]};
  const auto je{task()._y_range[1]};
  const auto ke{task()._z_range[1]};

  auto emf{emfPtr()};
  const auto ca_x{cax()};
  const auto ca_y{cay()};
  const auto ca_z{caz()};

  // xn
  const auto js_xn = jStartXN();
  const auto je_xn = jEndXN();
  const auto ks_xn = kStartXN();
  const auto ke_xn = kEndXN();
  for (auto j{js_xn}; j < je_xn; ++j) {
    for (auto k{ks_xn}; k < ke_xn + 1; ++k) {
      auto& ey_xn{emfPtr()->ey()(is, j, k)};
      const auto hz_i{hzInc(is - 1, j, k)};
      ey_xn += ca_x * hz_i;
    }
  }
  for (auto j{js_xn}; j < je_xn + 1; ++j) {
    for (auto k{ks_xn}; k < ke_xn; ++k) {
      auto& ez_xn{emfPtr()->ez()(is, j, k)};
      const auto hy_i{hyInc(is - 1, j, k)};
      ez_xn -= ca_x * hy_i;
    }
  }
  // xp
  const auto js_xp = jStartXP();
  const auto je_xp = jEndXP();
  const auto ks_xp = kStartXP();
  const auto ke_xp = kEndXP();
  for (auto j{js_xp}; j < je_xp; ++j) {
    for (auto k{ks_xp}; k < ke_xp + 1; ++k) {
      auto& ey_xp{emfPtr()->ey()(ie, j, k)};
      const auto hz_i{hzInc(ie, j, k)};
      ey_xp -= ca_x * hz_i;
    }
  }
  for (auto j{js_xp}; j < je_xp + 1; ++j) {
    for (auto k{ks_xp}; k < ke_xp; ++k) {
      auto& ez_xp{emfPtr()->ez()(ie, j, k)};
      const auto hy_i{hyInc(ie, j, k)};
      ez_xp += ca_x * hy_i;
    }
  }

  // yn
  const auto is_yn = iStartYN();
  const auto ie_yn = iEndYN();
  const auto ks_yn = kStartYN();
  const auto ke_yn = kEndYN();
  for (auto i{is_yn}; i < ie_yn + 1; ++i) {
    for (auto k{ks_yn}; k < ke_yn; ++k) {
      auto& ez_yn{emfPtr()->ez()(i, js, k)};
      const auto hx_i{hxInc(i, js - 1, k)};
      ez_yn += ca_y * hx_i;
    }
  }
  for (auto i{is_yn}; i < ie_yn; ++i) {
    for (auto k{ks_yn}; k < ke_yn + 1; ++k) {
      auto& ex_yn{emfPtr()->ex()(i, js, k)};
      const auto hz_i{hzInc(i, js - 1, k)};
      ex_yn -= ca_y * hz_i;
    }
  }
  // yp
  const auto is_yp = iStartYP();
  const auto ie_yp = iEndYP();
  const auto ks_yp = kStartYP();
  const auto ke_yp = kEndYP();
  for (auto i{is_yp}; i < ie_yp + 1; ++i) {
    for (auto k{ks_yp}; k < ke_yp; ++k) {
      auto& ez_yp{emfPtr()->ez()(i, je, k)};
      const auto hx_i{hxInc(i, je, k)};
      ez_yp -= ca_y * hx_i;
    }
  }
  for (auto i{is_yp}; i < ie_yp; ++i) {
    for (auto k{ks_yp}; k < ke_yp + 1; ++k) {
      auto& ex_yp{emfPtr()->ex()(i, je, k)};
      const auto hz_i{hzInc(i, je, k)};
      ex_yp += ca_y * hz_i;
    }
  }

  // zn
  const auto is_zn = iStartZN();
  const auto ie_zn = iEndZN();
  const auto js_zn = jStartZN();
  const auto je_zn = jEndZN();
  for (auto i{is_zn}; i < ie_zn; ++i) {
    for (auto j{js_zn}; j < je_zn + 1; ++j) {
      auto& ex_zn{emfPtr()->ex()(i, j, ks)};
      const auto hy_i{hyInc(i, j, ks - 1)};
      ex_zn += ca_z * hy_i;
    }
  }
  for (auto i{is_zn}; i < ie_zn + 1; ++i) {
    for (auto j{js_zn}; j < je_zn; ++j) {
      auto& ey_zn{emfPtr()->ey()(i, j, ks)};
      const auto hx_i{hxInc(i, j, ks - 1)};
      ey_zn -= ca_z * hx_i;
    }
  }
  // zp
  const auto is_zp = iStartZP();
  const auto ie_zp = iEndZP();
  const auto js_zp = jStartZP();
  const auto je_zp = jEndZP();
  for (auto i{is_zp}; i < ie_zp; ++i) {
    for (auto j{js_zp}; j < je_zp + 1; ++j) {
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

void TFSF3DCorrector::correctH() {
  const auto is{task()._x_range[0]};
  const auto js{task()._y_range[0]};
  const auto ks{task()._z_range[0]};
  const auto ie{task()._x_range[1]};
  const auto je{task()._y_range[1]};
  const auto ke{task()._z_range[1]};

  auto emf{emfPtr()};
  const auto cb_x{cbx()};
  const auto cb_y{cby()};
  const auto cb_z{cbz()};

  // xn
  const auto js_xn = jStartXN();
  const auto je_xn = jEndXN();
  const auto ks_xn = kStartXN();
  const auto ke_xn = kEndXN();
  for (auto j{js_xn}; j < je_xn; ++j) {
    for (auto k{ks_xn}; k < ke_xn + 1; ++k) {
      auto& hz_xn{emfPtr()->hz()(is - 1, j, k)};
      const auto ey_i{eyInc(is, j, k)};
      hz_xn += cb_x * ey_i;
    }
  }
  for (auto j{js_xn}; j < je_xn + 1; ++j) {
    for (auto k{ks_xn}; k < ke_xn; ++k) {
      auto& hy_xn{emfPtr()->hy()(is - 1, j, k)};
      const auto ez_i{ezInc(is, j, k)};
      hy_xn -= cb_x * ez_i;
    }
  }
  // xp
  const auto js_xp = jStartXP();
  const auto je_xp = jEndXP();
  const auto ks_xp = kStartXP();
  const auto ke_xp = kEndXP();
  for (auto j{js_xp}; j < je_xp; ++j) {
    for (auto k{ks_xp}; k < ke_xp + 1; ++k) {
      auto& hz_xp{emfPtr()->hz()(ie, j, k)};
      const auto ey_i{eyInc(ie, j, k)};
      hz_xp -= cb_x * ey_i;
    }
  }
  for (auto j{js_xp}; j < je_xp + 1; ++j) {
    for (auto k{ks_xp}; k < ke_xp; ++k) {
      auto& hy_xp{emfPtr()->hy()(ie, j, k)};
      const auto ez_i{ezInc(ie, j, k)};
      hy_xp += cb_x * ez_i;
    }
  }

  // yn
  const auto is_yn = iStartYN();
  const auto ie_yn = iEndYN();
  const auto ks_yn = kStartYN();
  const auto ke_yn = kEndYN();
  for (auto i{is_yn}; i < ie_yn + 1; ++i) {
    for (auto k{ks_yn}; k < ke_yn; ++k) {
      auto& hx_yn{emfPtr()->hx()(i, js - 1, k)};
      const auto ez_i{ezInc(i, js, k)};
      hx_yn += cb_y * ez_i;
    }
  }
  for (auto i{is_yn}; i < ie_yn; ++i) {
    for (auto k{ks_yn}; k < ke_yn + 1; ++k) {
      auto& hz_yn{emfPtr()->hz()(i, js - 1, k)};
      const auto ex_i{exInc(i, js, k)};
      hz_yn -= cb_y * ex_i;
    }
  }
  // yp
  const auto is_yp = iStartYP();
  const auto ie_yp = iEndYP();
  const auto ks_yp = kStartYP();
  const auto ke_yp = kEndYP();
  for (auto i{is_yp}; i < ie_yp + 1; ++i) {
    for (auto k{ks_yp}; k < ke_yp; ++k) {
      auto& hx_yp{emfPtr()->hx()(i, je, k)};
      const auto ez_i{ezInc(i, je, k)};
      hx_yp -= cb_y * ez_i;
    }
  }
  for (auto i{is_yp}; i < ie_yp; ++i) {
    for (auto k{ks_yp}; k < ke_yp + 1; ++k) {
      auto& hz_yp{emfPtr()->hz()(i, je, k)};
      const auto ex_i{exInc(i, je, k)};
      hz_yp += cb_y * ex_i;
    }
  }

  // zn
  const auto is_zn = iStartZN();
  const auto ie_zn = iEndZN();
  const auto js_zn = jStartZN();
  const auto je_zn = jEndZN();
  for (auto i{is_zn}; i < ie_zn; ++i) {
    for (auto j{js_zn}; j < je_zn + 1; ++j) {
      auto& hy_zn{emfPtr()->hy()(i, j, ks - 1)};
      const auto ex_i{exInc(i, j, ks)};
      hy_zn += cb_z * ex_i;
    }
  }
  for (auto i{is_zn}; i < ie_zn + 1; ++i) {
    for (auto j{js_zn}; j < je_zn; ++j) {
      auto& hx_zn{emfPtr()->hx()(i, j, ks - 1)};
      const auto ey_i{eyInc(i, j, ks)};
      hx_zn -= cb_z * ey_i;
    }
  }
  // zp
  const auto is_zp = iStartZP();
  const auto ie_zp = iEndZP();
  const auto js_zp = jStartZP();
  const auto je_zp = jEndZP();
  for (auto i{is_zp}; i < ie_zp; ++i) {
    for (auto j{js_zp}; j < je_zp + 1; ++j) {
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

}  // namespace xfdtd
