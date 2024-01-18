#ifndef _XFDTD_LIB_ELECTROMAGNETIC_FIELD_H_
#define _XFDTD_LIB_ELECTROMAGNETIC_FIELD_H_

#include <xfdtd/exception/exception.h>

#include <utility>
#include <xtensor/xarray.hpp>

#include "xfdtd/coordinate_system/coordinate_system.h"

namespace xfdtd {

class XFDTDEMFException : public XFDTDException {
 public:
  explicit XFDTDEMFException(std::string message = "XFDTD EMF Exception")
      : XFDTDException(std::move(message)) {}
};

class EMF {
 public:
  enum class Component { X, Y, Z, Magnitude };

  enum class Attribute { E, H, J };

  enum class Field { EX, EY, EZ, EM, HX, HY, HZ, HM, JX, JY, JZ, JM };

  static Field fieldFromAttributeAndComponent(Attribute a, Component c);

  static Component componentFromField(Field f);

  static Attribute attributeFromField(Field f);

  EMF() = default;

  EMF(const EMF&) = default;

  EMF(EMF&&) noexcept = default;

  EMF& operator=(const EMF&) = default;

  EMF& operator=(EMF&&) noexcept = default;

  ~EMF() = default;

  const xt::xarray<double>& ex() const;

  const xt::xarray<double>& ey() const;

  const xt::xarray<double>& ez() const;

  const xt::xarray<double>& hx() const;

  const xt::xarray<double>& hy() const;

  const xt::xarray<double>& hz() const;

  const xt::xarray<double>& jx() const;

  const xt::xarray<double>& jy() const;

  const xt::xarray<double>& jz() const;

  const xt::xarray<double>& jxPrev() const;

  const xt::xarray<double>& jyPrev() const;

  const xt::xarray<double>& jzPrev() const;

  const xt::xarray<double>& exPrev() const;

  const xt::xarray<double>& eyPrev() const;

  const xt::xarray<double>& ezPrev() const;

  const xt::xarray<double>& field(Field f) const;

  const xt::xarray<double>& field(Attribute a, Component c) const;

  xt::xarray<double>& ex();

  xt::xarray<double>& ey();

  xt::xarray<double>& ez();

  xt::xarray<double>& hx();

  xt::xarray<double>& hy();

  xt::xarray<double>& hz();

  xt::xarray<double>& jx();

  xt::xarray<double>& jy();

  xt::xarray<double>& jz();

  xt::xarray<double>& jxPrev();

  xt::xarray<double>& jyPrev();

  xt::xarray<double>& jzPrev();

  xt::xarray<double>& exPrev();

  xt::xarray<double>& eyPrev();

  xt::xarray<double>& ezPrev();

  xt::xarray<double>& field(Field f);

  xt::xarray<double>& field(Attribute a, Component c);

  double fieldFaceCenter(std::size_t i, std::size_t j, std::size_t k, Field f,
                         Axis::XYZ xyz) const;

  double fieldFaceCenterX(std::size_t i, std::size_t j, std::size_t k,
                          Field f) const;

  double fieldFaceCenterY(std::size_t i, std::size_t j, std::size_t k,
                          Field f) const;

  double fieldFaceCenterZ(std::size_t i, std::size_t j, std::size_t k,
                          Field f) const;

  // double exFaceCenterX(std::size_t i, std::size_t j, std::size_t k) const;

  double eyFaceCenterX(std::size_t i, std::size_t j, std::size_t k) const;

  double ezFaceCenterX(std::size_t i, std::size_t j, std::size_t k) const;

  // double hxFaceCenterX(std::size_t i, std::size_t j, std::size_t k) const;

  double hyFaceCenterX(std::size_t i, std::size_t j, std::size_t k) const;

  double hzFaceCenterX(std::size_t i, std::size_t j, std::size_t k) const;

  double ezFaceCenterY(std::size_t i, std::size_t j, std::size_t k) const;

  double exFaceCenterY(std::size_t i, std::size_t j, std::size_t k) const;

  double hzFaceCenterY(std::size_t i, std::size_t j, std::size_t k) const;

  double hxFaceCenterY(std::size_t i, std::size_t j, std::size_t k) const;

  double exFaceCenterZ(std::size_t i, std::size_t j, std::size_t k) const;

  double eyFaceCenterZ(std::size_t i, std::size_t j, std::size_t k) const;

  double hxFaceCenterZ(std::size_t i, std::size_t j, std::size_t k) const;

  double hyFaceCenterZ(std::size_t i, std::size_t j, std::size_t k) const;

  void allocateEx(std::size_t nx, std::size_t ny, std::size_t nz);

  void allocateEy(std::size_t nx, std::size_t ny, std::size_t nz);

  void allocateEz(std::size_t nx, std::size_t ny, std::size_t nz);

  void allocateHx(std::size_t nx, std::size_t ny, std::size_t nz);

  void allocateHy(std::size_t nx, std::size_t ny, std::size_t nz);

  void allocateHz(std::size_t nx, std::size_t ny, std::size_t nz);

  void allocateJx(std::size_t nx, std::size_t ny, std::size_t nz);

  void allocateJy(std::size_t nx, std::size_t ny, std::size_t nz);

  void allocateJz(std::size_t nx, std::size_t ny, std::size_t nz);

  void allocateJxPrev(std::size_t nx, std::size_t ny, std::size_t nz);

  void allocateJyPrev(std::size_t nx, std::size_t ny, std::size_t nz);

  void allocateJzPrev(std::size_t nx, std::size_t ny, std::size_t nz);

  void allocateExPrev(std::size_t nx, std::size_t ny, std::size_t nz);

  void allocateEyPrev(std::size_t nx, std::size_t ny, std::size_t nz);

  void allocateEzPrev(std::size_t nx, std::size_t ny, std::size_t nz);

 private:
  xt::xarray<double> _ex, _ey, _ez;
  xt::xarray<double> _hx, _hy, _hz;
  // polarization current
  xt::xarray<double> _jx, _jy, _jz;
  // polarization current arrays at previous time step
  xt::xarray<double> _jx_prev, _jy_prev, _jz_prev;
  // E-field (at previous time step) for Lorentz only.
  xt::xarray<double> _ex_prev, _ey_prev, _ez_prev;
};

}  // namespace xfdtd

#endif  // _XFDTD_LIB_ELECTROMAGNETIC_FIELD_H_
