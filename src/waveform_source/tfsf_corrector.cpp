#include "waveform_source/tfsf_corrector.h"

#include <memory>
#include <sstream>

namespace xfdtd {

std::string TFSFCorrector::toString() const {
  std::stringstream ss;
  ss << "Total Field Scattered Field Corrector:\n";
  ss << " XN: "
     << "EY: " << globalEyTaskXN().toString()
     << " EZ: " << globalEzTaskXN().toString() << "\n";
  ss << " XP: "
     << "EY: " << globalEyTaskXP().toString()
     << " EZ: " << globalEzTaskXP().toString() << "\n";
  ss << " YN: "
     << "EZ: " << globalEzTaskYN().toString()
     << " EX: " << globalExTaskYN().toString() << "\n";
  ss << " YP: "
     << "EZ: " << globalEzTaskYP().toString()
     << " EX: " << globalExTaskYP().toString() << "\n";
  ss << " ZN: "
     << "EX: " << globalExTaskZN().toString()
     << " EY: " << globalEyTaskZN().toString() << "\n";
  ss << " ZP: "
     << "EX: " << globalExTaskZP().toString()
     << " EY: " << globalEyTaskZP().toString();
  return ss.str();
}

Real TFSFCorrector::exInc(std::size_t i, std::size_t j, std::size_t k) const {
  i = i - globalStart().i() + 1;
  j = j - globalStart().j();
  k = k - globalStart().k();
  auto projection{_projection_x_half(i) + _projection_y_int(j) +
                  _projection_z_int(k)};
  auto index{static_cast<std::size_t>(projection)};
  auto weight{projection - index};
  return (1 - weight) * _ex_inc(index) + weight * _ex_inc(index + 1);
}

Real TFSFCorrector::eyInc(std::size_t i, std::size_t j, std::size_t k) const {
  i = i - globalStart().i();
  j = j - globalStart().j() + 1;
  k = k - globalStart().k();
  auto projection{_projection_x_int(i) + _projection_y_half(j) +
                  _projection_z_int(k)};
  auto index{static_cast<std::size_t>(projection)};
  auto weight{projection - index};
  return (1 - weight) * _ey_inc(index) + weight * _ey_inc(index + 1);
}

Real TFSFCorrector::ezInc(std::size_t i, std::size_t j, std::size_t k) const {
  i = i - globalStart().i();
  j = j - globalStart().j();
  k = k - globalStart().k() + 1;
  auto projection{_projection_x_int(i) + _projection_y_int(j) +
                  _projection_z_half(k)};
  auto index{static_cast<std::size_t>(projection)};
  auto weight{projection - index};
  return (1 - weight) * _ez_inc(index) + weight * _ez_inc(index + 1);
}

Real TFSFCorrector::hxInc(std::size_t i, std::size_t j, std::size_t k) const {
  i = i - globalStart().i();
  j = j - globalStart().j() + 1;
  k = k - globalStart().k() + 1;
  auto projection{_projection_x_int(i) + _projection_y_half(j) +
                  _projection_z_half(k) - 0.5};
  auto index{static_cast<std::size_t>(projection)};
  auto weight{projection - index};
  return (1 - weight) * _hx_inc(index) + weight * _hx_inc(index + 1);
}

Real TFSFCorrector::hyInc(std::size_t i, std::size_t j, std::size_t k) const {
  i = i - globalStart().i() + 1;
  j = j - globalStart().j();
  k = k - globalStart().k() + 1;
  auto projection{_projection_x_half(i) + _projection_y_int(j) +
                  _projection_z_half(k) - 0.5};
  auto index{static_cast<std::size_t>(projection)};
  auto weight{projection - index};
  return (1 - weight) * _hy_inc(index) + weight * _hy_inc(index + 1);
}

Real TFSFCorrector::hzInc(std::size_t i, std::size_t j, std::size_t k) const {
  i = i - globalStart().i() + 1;
  j = j - globalStart().j() + 1;
  k = k - globalStart().k();
  auto projection{_projection_x_half(i) + _projection_y_half(j) +
                  _projection_z_int(k) - 0.5};
  auto index{static_cast<std::size_t>(projection)};
  auto weight{projection - index};
  return (1 - weight) * _hz_inc(index) + weight * _hz_inc(index + 1);
}

void TFSFCorrector::correctEyXN() {
  const auto& task = globalEyTaskXN();
  const auto& offset = gridSpace()->globalBox().origin();

  const auto is = task.xRange().start();
  const auto js = task.yRange().start();
  const auto je = task.yRange().end();
  const auto ks = task.zRange().start();
  const auto ke = task.zRange().end();

  auto&& emf = emfPtr();
  const auto ca_x{cax()};

  const auto is_node = is - offset.i();
  for (auto j{js}; j < js; ++j) {
    const auto j_node = j - offset.j();
    for (auto k{ks}; k < ke; ++k) {
      const auto k_node = k - offset.k();
      auto& ey_xn{emf->ey()(is_node, j_node, k_node)};
      const auto hz_i{hzInc(is - 1, j, k)};
      ey_xn += ca_x * hz_i;
    }
  }
}

void TFSFCorrector::correctEzXN() {
  const auto& task = globalEzTaskXN();
  const auto& offset = gridSpace()->globalBox().origin();

  const auto is = task.xRange().start();
  const auto js = task.yRange().start();
  const auto je = task.yRange().end();
  const auto ks = task.zRange().start();
  const auto ke = task.zRange().end();

  auto&& emf = emfPtr();
  const auto ca_x{cax()};

  const auto is_node = is - offset.i();
  for (auto j{js}; j < je; ++j) {
    const auto j_node = j - offset.j();
    for (auto k{ks}; k < ke; ++k) {
      const auto k_node = k - offset.k();
      auto& ez_xn{emf->ez()(is_node, j_node, k_node)};
      const auto hy_i{hyInc(is - 1, j, k)};
      ez_xn -= ca_x * hy_i;
    }
  }
}

void TFSFCorrector::correctEyXP() {
  const auto& task = globalEyTaskXP();
  const auto& offset = gridSpace()->globalBox().origin();

  const auto ie = task.xRange().start();  // careful!!
  const auto js = task.yRange().start();
  const auto je = task.yRange().end();
  const auto ks = task.zRange().start();
  const auto ke = task.zRange().end();

  auto&& emf = emfPtr();
  const auto ca_x{cax()};

  const auto ie_node = ie - offset.i();
  for (auto j{js}; j < je; ++j) {
    const auto j_node = j - offset.j();
    for (auto k{ks}; k < ke; ++k) {
      const auto k_node = k - offset.k();
      auto& ey_xp{emf->ey()(ie_node, j_node, k_node)};
      const auto hz_i{hzInc(ie, j, k)};
      ey_xp -= ca_x * hz_i;
    }
  }
}

void TFSFCorrector::correctEzXP() {
  const auto& task = globalEzTaskXP();
  const auto& offset = gridSpace()->globalBox().origin();

  const auto ie = task.xRange().start();
  const auto js = task.yRange().start();
  const auto je = task.yRange().end();
  const auto ks = task.zRange().start();
  const auto ke = task.zRange().end();

  auto&& emf = emfPtr();
  const auto ca_x{cax()};

  const auto ie_node = ie - offset.i();
  for (auto j{js}; j < je; ++j) {
    const auto j_node = j - offset.j();
    for (auto k{ks}; k < ke; ++k) {
      const auto k_node = k - offset.k();
      auto& ez_xp{emf->ez()(ie_node, j_node, k_node)};
      const auto hy_i{hyInc(ie, j, k)};
      ez_xp += ca_x * hy_i;
    }
  }
}

void TFSFCorrector::correctEzYN() {
  const auto& task = globalEzTaskYN();
  const auto& offset = gridSpace()->globalBox().origin();

  const auto is = task.xRange().start();
  const auto ie = task.xRange().end();
  const auto js = task.yRange().start();
  const auto ks = task.zRange().start();
  const auto ke = task.zRange().end();

  auto&& emf = emfPtr();
  const auto ca_y = cay();

  const auto js_node = js - offset.j();
  for (auto i{is}; i < ie; ++i) {
    const auto i_node = i - offset.i();
    for (auto k{ks}; k < ke; ++k) {
      const auto k_node = k - offset.k();
      auto& ez_yn{emf->ez()(i_node, js_node, k_node)};
      const auto hx_i{hxInc(i, js - 1, k)};
      ez_yn += ca_y * hx_i;
    }
  }
}

void TFSFCorrector::correctExYN() {
  const auto& task = globalExTaskYN();
  const auto& offset = gridSpace()->globalBox().origin();

  const auto is = task.xRange().start();
  const auto ie = task.xRange().end();
  const auto js = task.yRange().start();
  const auto ks = task.zRange().start();
  const auto ke = task.zRange().end();

  auto&& emf = emfPtr();
  const auto ca_y = cay();

  const auto js_node = js - offset.j();
  for (auto i{is}; i < ie; ++i) {
    const auto i_node = i - offset.i();
    for (auto k{ks}; k < ke; ++k) {
      const auto k_node = k - offset.k();
      auto& ex_yn{emf->ex()(i_node, js_node, k_node)};
      const auto hz_i{hzInc(i, js - 1, k)};
      ex_yn -= ca_y * hz_i;
    }
  }
}

void TFSFCorrector::correctEzYP() {
  const auto& task = globalEzTaskYP();
  const auto& offset = gridSpace()->globalBox().origin();

  const auto is = task.xRange().start();
  const auto ie = task.xRange().end();
  const auto je = task.yRange().start();
  const auto ks = task.zRange().start();
  const auto ke = task.zRange().end();

  auto&& emf = emfPtr();
  const auto ca_y = cay();

  const auto je_node = je - offset.j();
  for (auto i{is}; i < ie; ++i) {
    const auto i_node = i - offset.i();
    for (auto k{ks}; k < ke; ++k) {
      const auto k_node = k - offset.k();
      auto& ez_yp{emf->ez()(i_node, je_node, k_node)};
      const auto hx_i{hxInc(i, je, k)};
      ez_yp -= ca_y * hx_i;
    }
  }
}

void TFSFCorrector::correctExYP() {
  const auto& task = globalExTaskYP();
  const auto& offset = gridSpace()->globalBox().origin();

  const auto is = task.xRange().start();
  const auto ie = task.xRange().end();
  const auto je = task.yRange().start();
  const auto ks = task.zRange().start();
  const auto ke = task.zRange().end();

  auto&& emf = emfPtr();
  const auto ca_y = cay();

  const auto je_node = je - offset.j();
  for (auto i{is}; i < ie; ++i) {
    const auto i_node = i - offset.i();
    for (auto k{ks}; k < ke; ++k) {
      const auto k_node = k - offset.k();
      auto& ex_yp{emf->ez()(i_node, je_node, k_node)};
      const auto hx_i{hxInc(i, je, k)};
      ex_yp += ca_y * hx_i;
    }
  }
}

void TFSFCorrector::correctExZN() {
  const auto& task = globalExTaskZN();
  const auto& offset = gridSpace()->globalBox().origin();

  const auto is = task.xRange().start();
  const auto ie = task.xRange().end();
  const auto js = task.yRange().start();
  const auto je = task.yRange().end();
  const auto ks = task.zRange().start();

  auto&& emf = emfPtr();
  const auto ca_z = caz();

  const auto ks_node = ks - offset.k();
  for (auto i{is}; i < ie; ++i) {
    const auto i_node = i - offset.i();
    for (auto j{js}; j < je; ++j) {
      const auto j_node = j - offset.j();
      auto& ex_zn{emf->ex()(i_node, j_node, ks_node)};
      const auto hy_i{hyInc(i, j, ks - 1)};
      ex_zn += ca_z * hy_i;
    }
  }
}

void TFSFCorrector::correctEyZN() {
  const auto& task = globalEyTaskZN();
  const auto& offset = gridSpace()->globalBox().origin();

  const auto is = task.xRange().start();
  const auto ie = task.xRange().end();
  const auto js = task.yRange().start();
  const auto je = task.yRange().end();
  const auto ks = task.zRange().start();

  auto&& emf = emfPtr();
  const auto ca_z = caz();

  const auto ks_node = ks - offset.k();
  for (auto i{is}; i < ie; ++i) {
    const auto i_node = i - offset.i();
    for (auto j{js}; j < je; ++j) {
      const auto j_node = j - offset.j();
      auto& ey_zn{emf->ey()(i_node, j_node, ks_node)};
      const auto hx_i{hxInc(i, j, ks - 1)};
      ey_zn -= ca_z * hx_i;
    }
  }
}

void TFSFCorrector::correctExZP() {
  const auto& task = globalExTaskZP();
  const auto& offset = gridSpace()->globalBox().origin();

  const auto is = task.xRange().start();
  const auto ie = task.xRange().end();
  const auto js = task.yRange().start();
  const auto je = task.yRange().end();
  const auto ke = task.zRange().start();

  auto&& emf = emfPtr();
  const auto ca_z = caz();

  const auto ke_node = ke - offset.k();
  for (auto i{is}; i < ie; ++i) {
    const auto i_node = i - offset.i();
    for (auto j{js}; j < je; ++j) {
      const auto j_node = j - offset.j();
      auto& ex_zp{emf->ex()(i_node, j_node, ke_node)};
      const auto hy_i{hyInc(i, j, ke)};
      ex_zp -= ca_z * hy_i;
    }
  }
}

void TFSFCorrector::correctEyZP() {
  const auto& task = globalEyTaskZP();
  const auto& offset = gridSpace()->globalBox().origin();

  const auto is = task.xRange().start();
  const auto ie = task.xRange().end();
  const auto js = task.yRange().start();
  const auto je = task.yRange().end();
  const auto ke = task.zRange().start();

  auto&& emf = emfPtr();
  const auto ca_z = caz();

  const auto ke_node = ke - offset.k();
  for (auto i{is}; i < ie; ++i) {
    const auto i_node = i - offset.i();
    for (auto j{js}; j < je; ++j) {
      const auto j_node = j - offset.j();
      auto& ey_zp{emf->ey()(i_node, j_node, ke_node)};
      const auto hx_i{hxInc(i, j, ke)};
      ey_zp += ca_z * hx_i;
    }
  }
}

void TFSFCorrector::correctHzXN() {
  const auto& task = globalEyTaskXN();
  const auto& offset = gridSpace()->globalBox().origin();

  const auto is = task.xRange().start();
  const auto js = task.yRange().start();
  const auto je = task.yRange().end();
  const auto ks = task.zRange().start();
  const auto ke = task.zRange().end();

  auto&& emf = emfPtr();
  const auto cb_x{cbx()};

  const auto is_node = is - offset.i();
  for (auto j{js}; j < je; ++j) {
    const auto j_node = j - offset.j();
    for (auto k{ks}; k < ke; ++k) {
      const auto k_node = k - offset.k();
      auto& hz_xn{emf->hz()(is_node - 1, j_node, k_node)};
      const auto ey_i{eyInc(is, j, k)};
      hz_xn += cb_x * ey_i;
    }
  }
}

void TFSFCorrector::correctHyXN() {
  const auto& task = globalEzTaskXN();
  const auto& offset = gridSpace()->globalBox().origin();

  const auto is = task.xRange().start();
  const auto js = task.yRange().start();
  const auto je = task.yRange().end();
  const auto ks = task.zRange().start();
  const auto ke = task.zRange().end();

  auto&& emf = emfPtr();
  const auto cb_x{cbx()};
  const auto is_node = is - offset.i();
  for (auto j{js}; j < je; ++j) {
    const auto j_node = j - offset.j();
    for (auto k{ks}; k < ke; ++k) {
      const auto k_node = k - offset.k();
      auto& hy_xn{emf->hy()(is_node - 1, j_node, k_node)};
      const auto ez_i{ezInc(is, j, k)};
      hy_xn -= cb_x * ez_i;
    }
  }
}

void TFSFCorrector::correctHzXP() {
  const auto& task = globalEyTaskXP();
  const auto& offset = gridSpace()->globalBox().origin();

  const auto ie = task.xRange().start();
  const auto js = task.yRange().start();
  const auto je = task.yRange().end();
  const auto ks = task.zRange().start();
  const auto ke = task.zRange().end();

  auto&& emf = emfPtr();
  const auto cb_x{cbx()};

  const auto ie_node = ie - offset.i();
  for (auto j{js}; j < je; ++j) {
    const auto j_node = j - offset.j();
    for (auto k{ks}; k < ke; ++k) {
      const auto k_node = k - offset.k();
      auto& hz_xp{emf->hz()(ie_node, j_node, k_node)};
      const auto ey_i{eyInc(ie, j, k)};
      hz_xp -= cb_x * ey_i;
    }
  }
}

void TFSFCorrector::correctHyXP() {
  const auto& task = globalEzTaskXP();
  const auto& offset = gridSpace()->globalBox().origin();

  const auto ie = task.xRange().start();
  const auto js = task.yRange().start();
  const auto je = task.yRange().end();
  const auto ks = task.zRange().start();
  const auto ke = task.zRange().end();

  auto&& emf = emfPtr();
  const auto cb_x{cbx()};

  const auto ie_node = ie - offset.i();
  for (auto j{js}; j < je; ++j) {
    const auto j_node = j - offset.j();
    for (auto k{ks}; k < ke; ++k) {
      const auto k_node = k - offset.k();
      auto& hy_xp{emf->hy()(ie_node, j_node, k_node)};
      const auto ez_i{ezInc(ie, j, k)};
      hy_xp += cb_x * ez_i;
    }
  }
}

void TFSFCorrector::correctHxYN() {
  const auto& task = globalEzTaskYN();
  const auto& offset = gridSpace()->globalBox().origin();

  const auto is = task.xRange().start();
  const auto ie = task.xRange().end();
  const auto js = task.yRange().start();
  const auto ks = task.zRange().start();
  const auto ke = task.zRange().end();

  auto&& emf = emfPtr();
  const auto cb_y{cby()};

  const auto js_node = js - offset.j();
  for (auto i{is}; i < ie; ++i) {
    const auto i_node = i - offset.i();
    for (auto k{ks}; k < ke; ++k) {
      const auto k_node = k - offset.k();
      auto& hx_yn{emf->hx()(i_node, js_node - 1, k_node)};
      const auto ez_i{ezInc(i, js, k)};
      hx_yn += cb_y * ez_i;
    }
  }
}

void TFSFCorrector::correctHzYN() {
  const auto& task = globalExTaskYN();
  const auto& offset = gridSpace()->globalBox().origin();

  const auto is = task.xRange().start();
  const auto ie = task.xRange().end();
  const auto js = task.yRange().start();
  const auto ks = task.zRange().start();
  const auto ke = task.zRange().end();

  auto&& emf = emfPtr();
  const auto cb_y{cby()};

  const auto js_node = js - offset.j();
  for (auto i{is}; i < ie; ++i) {
    const auto i_node = i - offset.i();
    for (auto k{ks}; k < ke; ++k) {
      const auto k_node = k - offset.k();
      auto& hz_yn{emf->hz()(i_node, js_node - 1, k_node)};
      const auto ex_i{exInc(i, js, k)};
      hz_yn -= cb_y * ex_i;
    }
  }
}

void TFSFCorrector::correctHxYP() {
  const auto& task = globalEzTaskYP();
  const auto& offset = gridSpace()->globalBox().origin();

  const auto is = task.xRange().start();
  const auto ie = task.xRange().end();
  const auto je = task.yRange().start();
  const auto ks = task.zRange().start();
  const auto ke = task.zRange().end();

  auto&& emf = emfPtr();
  const auto cb_y{cby()};

  const auto je_node = je - offset.j();
  for (auto i{is}; i < ie; ++i) {
    const auto i_node = i - offset.i();
    for (auto k{ks}; k < ke; ++k) {
      const auto k_node = k - offset.k();
      auto& hx_yp{emf->hx()(i_node, je_node, k_node)};
      const auto ez_i{ezInc(i, je, k)};
      hx_yp -= cb_y * ez_i;
    }
  }
}

void TFSFCorrector::correctHzYP() {
  const auto& task = globalExTaskYP();
  const auto& offset = gridSpace()->globalBox().origin();

  const auto is = task.xRange().start();
  const auto ie = task.xRange().end();
  const auto je = task.yRange().start();
  const auto ks = task.zRange().start();
  const auto ke = task.zRange().end();

  auto&& emf = emfPtr();
  const auto cb_y{cby()};

  const auto je_node = je - offset.j();
  for (auto i{is}; i < ie; ++i) {
    const auto i_node = i - offset.i();
    for (auto k{ks}; k < ke; ++k) {
      const auto k_node = k - offset.k();
      auto& hz_yp{emf->hz()(i_node, je_node, k_node)};
      const auto ex_i{exInc(i, je, k)};
      hz_yp += cb_y * ex_i;
    }
  }
}

void TFSFCorrector::correctHxZN() {
  const auto& task = globalEyTaskZN();
  const auto& offset = gridSpace()->globalBox().origin();

  const auto is = task.xRange().start();
  const auto ie = task.xRange().end();
  const auto js = task.yRange().start();
  const auto je = task.yRange().end();
  const auto ks = task.zRange().start();

  auto&& emf = emfPtr();
  const auto cb_z{cbz()};

  const auto ks_node = ks - offset.k();
  for (auto i{is}; i < ie; ++i) {
    const auto i_node = i - offset.i();
    for (auto j{js}; j < je; ++j) {
      const auto j_node = j - offset.j();
      auto& hx_zn{emf->hx()(i_node, j_node, ks_node - 1)};
      const auto ey_i{eyInc(i, j, ks)};
      hx_zn -= cb_z * ey_i;
    }
  }
}

void TFSFCorrector::correctHyZN() {
  const auto& task = globalExTaskZN();
  const auto& offset = gridSpace()->globalBox().origin();

  const auto is = task.xRange().start();
  const auto ie = task.xRange().end();
  const auto js = task.yRange().start();
  const auto je = task.yRange().end();
  const auto ks = task.zRange().start();

  auto&& emf = emfPtr();
  const auto cb_z{cbz()};

  const auto ks_node = ks - offset.k();
  for (auto i{is}; i < ie; ++i) {
    const auto i_node = i - offset.i();
    for (auto j{js}; j < je; ++j) {
      const auto j_node = j - offset.j();
      auto& hy_zn{emf->hy()(i_node, j_node, ks_node - 1)};
      const auto ex_i{exInc(i, j, ks)};
      hy_zn += cb_z * ex_i;
    }
  }
}

void TFSFCorrector::correctHxZP() {
  const auto& task = globalEyTaskZP();
  const auto& offset = gridSpace()->globalBox().origin();

  const auto is = task.xRange().start();
  const auto ie = task.xRange().end();
  const auto js = task.yRange().start();
  const auto je = task.yRange().end();
  const auto ke = task.zRange().start();

  auto&& emf = emfPtr();
  const auto cb_z{cbz()};

  const auto ke_node = ke - offset.k();
  for (auto i{is}; i < ie; ++i) {
    const auto i_node = i - offset.i();
    for (auto j{js}; j < je; ++j) {
      const auto j_node = j - offset.j();
      auto& hx_zp{emf->hx()(i_node, j_node, ke_node)};
      const auto ey_i{eyInc(i, j, ke)};
      hx_zp += cb_z * ey_i;
    }
  }
}

void TFSFCorrector::correctHyZP() {
  const auto& task = globalExTaskZP();
  const auto& offset = gridSpace()->globalBox().origin();

  const auto is = task.xRange().start();
  const auto ie = task.xRange().end();
  const auto js = task.yRange().start();
  const auto je = task.yRange().end();
  const auto ke = task.zRange().start();

  auto&& emf = emfPtr();
  const auto cb_z{cbz()};

  const auto ke_node = ke - offset.k();
  for (auto i{is}; i < ie; ++i) {
    const auto i_node = i - offset.i();
    for (auto j{js}; j < je; ++j) {
      const auto j_node = j - offset.j();
      auto& hy_zp{emf->hy()(i_node, j_node, ke_node)};
      const auto ex_i{exInc(i, j, ke)};
      hy_zp -= cb_z * ex_i;
    }
  }
}

IndexTask TFSFCorrector::globalEyTaskXN() const {
  return _global_ey_task_xn;
}

IndexTask TFSFCorrector::globalEzTaskXN() const {
  return _global_ez_task_xn;
}

IndexTask TFSFCorrector::globalEyTaskXP() const {
  return _global_ey_task_xp;
}

IndexTask TFSFCorrector::globalEzTaskXP() const {
  return _global_ez_task_xp;
}

IndexTask TFSFCorrector::globalEzTaskYN() const {
  return _global_ez_task_yn;
}

IndexTask TFSFCorrector::globalExTaskYN() const {
  return _global_ex_task_yn;
}

IndexTask TFSFCorrector::globalEzTaskYP() const {
  return _global_ez_task_yp;
}

IndexTask TFSFCorrector::globalExTaskYP() const {
  return _global_ex_task_yp;
}

IndexTask TFSFCorrector::globalExTaskZN() const {
  return _global_ex_task_zn;
}

IndexTask TFSFCorrector::globalEyTaskZN() const {
  return _global_ey_task_zn;
}

IndexTask TFSFCorrector::globalExTaskZP() const {
  return _global_ex_task_zp;
}

IndexTask TFSFCorrector::globalEyTaskZP() const {
  return _global_ey_task_zp;
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
  correctEzXN();
  correctEzXP();
  correctEzYN();
  correctEzYP();
}

void TFSF2DCorrector::correctH() {
  correctHyXN();
  correctHyXP();
  correctHxYN();
  correctHxYP();
}

void TFSF3DCorrector::correctE() {
  correctEyXN();
  correctEzXN();
  correctEyXP();
  correctEzXP();
  correctEzYN();
  correctExYN();
  correctEzYP();
  correctExYP();
  correctExZN();
  correctEyZN();
  correctExZP();
  correctEyZP();
}

void TFSF3DCorrector::correctH() {
  correctHyXN();
  correctHzXN();
  correctHyXP();
  correctHzXP();
  correctHxYN();
  correctHzYN();
  correctHxYP();
  correctHzYP();
  correctHxZN();
  correctHyZN();
  correctHxZP();
  correctHyZP();
}

}  // namespace xfdtd
