#include "xfdtd/electromagnetic_field/electromagnetic_field.h"

#include <cstddef>

#include "xfdtd/coordinate_system/coordinate_system.h"

namespace xfdtd {

EMF::Field EMF::fieldFromAttributeAndComponent(EMF::Attribute a,
                                               EMF::Component c) {
  switch (a) {
    case EMF::Attribute::E:
      switch (c) {
        case EMF::Component::X:
          return EMF::Field::EX;
        case EMF::Component::Y:
          return EMF::Field::EY;
        case EMF::Component::Z:
          return EMF::Field::EZ;
        default:
          throw XFDTDEMFException("Invalid component type");
      }
    case EMF::Attribute::H:
      switch (c) {
        case EMF::Component::X:
          return EMF::Field::HX;
        case EMF::Component::Y:
          return EMF::Field::HY;
        case EMF::Component::Z:
          return EMF::Field::HZ;
        default:
          throw XFDTDEMFException("Invalid component type");
      }
    case EMF::Attribute::J:
      switch (c) {
        case EMF::Component::X:
          return EMF::Field::JX;
        case EMF::Component::Y:
          return EMF::Field::JY;
        case EMF::Component::Z:
          return EMF::Field::JZ;
        default:
          throw XFDTDEMFException("Invalid component type");
      }
    default:
      throw XFDTDEMFException("Invalid attribute type");
  }
}

EMF::Component EMF::componentFromField(EMF::Field f) {
  switch (f) {
    case EMF::Field::EX:
    case EMF::Field::HX:
    case EMF::Field::JX:
      return EMF::Component::X;
    case EMF::Field::EY:
    case EMF::Field::HY:
    case EMF::Field::JY:
      return EMF::Component::Y;
    case EMF::Field::EZ:
    case EMF::Field::HZ:
    case EMF::Field::JZ:
      return EMF::Component::Z;
    default:
      throw XFDTDEMFException("Invalid field type");
  }
}

EMF::Attribute EMF::attributeFromField(EMF::Field f) {
  switch (f) {
    case EMF::Field::EX:
    case EMF::Field::EY:
    case EMF::Field::EZ:
      return EMF::Attribute::E;
    case EMF::Field::HX:
    case EMF::Field::HY:
    case EMF::Field::HZ:
      return EMF::Attribute::H;
    case EMF::Field::JX:
    case EMF::Field::JY:
    case EMF::Field::JZ:
      return EMF::Attribute::J;
    default:
      throw XFDTDEMFException("Invalid field type");
  }
}

const xt::xarray<double>& EMF::ex() const { return _ex; }

const xt::xarray<double>& EMF::ey() const { return _ey; }

const xt::xarray<double>& EMF::ez() const { return _ez; }

const xt::xarray<double>& EMF::hx() const { return _hx; }

const xt::xarray<double>& EMF::hy() const { return _hy; }

const xt::xarray<double>& EMF::hz() const { return _hz; }

const xt::xarray<double>& EMF::jx() const { return _jx; }

const xt::xarray<double>& EMF::jy() const { return _jy; }

const xt::xarray<double>& EMF::jz() const { return _jz; }

const xt::xarray<double>& EMF::jxPrev() const { return _jx_prev; }

const xt::xarray<double>& EMF::jyPrev() const { return _jy_prev; }

const xt::xarray<double>& EMF::jzPrev() const { return _jz_prev; }

const xt::xarray<double>& EMF::exPrev() const { return _ex_prev; }

const xt::xarray<double>& EMF::eyPrev() const { return _ey_prev; }

const xt::xarray<double>& EMF::ezPrev() const { return _ez_prev; }

const xt::xarray<double>& EMF::field(Field f) const {
  switch (f) {
    case Field::EX:
      return _ex;
    case Field::EY:
      return _ey;
    case Field::EZ:
      return _ez;
    case Field::HX:
      return _hx;
    case Field::HY:
      return _hy;
    case Field::HZ:
      return _hz;
    case Field::JX:
      return _jx;
    case Field::JY:
      return _jy;
    case Field::JZ:
      return _jz;
    default:
      throw XFDTDEMFException("Invalid field type");
  }
}

const xt::xarray<double>& EMF::field(Attribute a, Component c) const {
  return field(EMF::fieldFromAttributeAndComponent(a, c));
}

xt::xarray<double>& EMF::ex() { return _ex; }

xt::xarray<double>& EMF::ey() { return _ey; }

xt::xarray<double>& EMF::ez() { return _ez; }

xt::xarray<double>& EMF::hx() { return _hx; }

xt::xarray<double>& EMF::hy() { return _hy; }

xt::xarray<double>& EMF::hz() { return _hz; }

xt::xarray<double>& EMF::jx() { return _jx; }

xt::xarray<double>& EMF::jy() { return _jy; }

xt::xarray<double>& EMF::jz() { return _jz; }

xt::xarray<double>& EMF::jxPrev() { return _jx_prev; }

xt::xarray<double>& EMF::jyPrev() { return _jy_prev; }

xt::xarray<double>& EMF::jzPrev() { return _jz_prev; }

xt::xarray<double>& EMF::exPrev() { return _ex_prev; }

xt::xarray<double>& EMF::eyPrev() { return _ey_prev; }

xt::xarray<double>& EMF::ezPrev() { return _ez_prev; }

xt::xarray<double>& EMF::field(Field f) {
  switch (f) {
    case Field::EX:
      return _ex;
    case Field::EY:
      return _ey;
    case Field::EZ:
      return _ez;
    case Field::HX:
      return _hx;
    case Field::HY:
      return _hy;
    case Field::HZ:
      return _hz;
    case Field::JX:
      return _jx;
    case Field::JY:
      return _jy;
    case Field::JZ:
      return _jz;
    default:
      throw XFDTDEMFException("Invalid field type");
  }
}

double EMF::fieldFaceCenter(std::size_t i, std::size_t j, std::size_t k,
                            Field f, Axis::XYZ xyz) const {
  switch (xyz) {
    case Axis::XYZ::X:
      return fieldFaceCenterX(i, j, k, f);
    case Axis::XYZ::Y:
      return fieldFaceCenterY(i, j, k, f);
    case Axis::XYZ::Z:
      return fieldFaceCenterZ(i, j, k, f);
    default:
      throw XFDTDEMFException("Invalid direction type");
  }
}

double EMF::fieldFaceCenterX(std::size_t i, std::size_t j, std::size_t k,
                             Field f) const {
  switch (f) {
    // case Field::EX:
    //   return exFaceCenterX(i, j, k);
    case Field::EY:
      return eyFaceCenterX(i, j, k);
    case Field::EZ:
      return ezFaceCenterX(i, j, k);
    // case Field::HX:
    //   return hxFaceCenterX(i, j, k);
    case Field::HY:
      return hyFaceCenterX(i, j, k);
    case Field::HZ:
      return hzFaceCenterX(i, j, k);
    default:
      throw XFDTDEMFException("Invalid field type");
  }
}

double EMF::fieldFaceCenterY(std::size_t i, std::size_t j, std::size_t k,
                             Field f) const {
  switch (f) {
    case Field::EX:
      return exFaceCenterY(i, j, k);
    case Field::EZ:
      return ezFaceCenterY(i, j, k);
    case Field::HX:
      return hxFaceCenterY(i, j, k);
    case Field::HZ:
      return hzFaceCenterY(i, j, k);
    default:
      throw XFDTDEMFException("Invalid field type");
  }
}

double EMF::fieldFaceCenterZ(std::size_t i, std::size_t j, std::size_t k,
                             Field f) const {
  switch (f) {
    case Field::EX:
      return exFaceCenterZ(i, j, k);
    case Field::EY:
      return eyFaceCenterZ(i, j, k);
    case Field::HX:
      return hxFaceCenterZ(i, j, k);
    case Field::HY:
      return hyFaceCenterZ(i, j, k);
    default:
      throw XFDTDEMFException("Invalid field type");
  }
}

double EMF::eyFaceCenterX(std::size_t i, std::size_t j, std::size_t k) const {
  return 0.5 * (_ey(i, j, k + 1) + _ey(i, j, k));
}

double EMF::ezFaceCenterX(std::size_t i, std::size_t j, std::size_t k) const {
  return 0.5 * (_ez(i, j + 1, k) + _ez(i, j, k));
}

// double EMF::hxFaceCenterX(std::size_t i, std::size_t j, std::size_t k) const
// {
//   return _hx(i, j, k);
// }

double EMF::hyFaceCenterX(std::size_t i, std::size_t j, std::size_t k) const {
  return 0.25 * (_hy(i, j + 1, k) + _hy(i, j, k) + _hy(i - 1, j + 1, k) +
                 _hy(i - 1, j, k));
}

double EMF::hzFaceCenterX(std::size_t i, std::size_t j, std::size_t k) const {
  return 0.25 * (_hz(i, j, k + 1) + _hz(i, j, k) + _hz(i - 1, j, k + 1) +
                 _hz(i - 1, j, k));
}

double EMF::ezFaceCenterY(std::size_t i, std::size_t j, std::size_t k) const {
  return 0.5 * (_ez(i + 1, j, k) + _ez(i, j, k));
}

double EMF::exFaceCenterY(std::size_t i, std::size_t j, std::size_t k) const {
  return 0.5 * (_ex(i, j, k + 1) + _ex(i, j, k));
}

double EMF::hzFaceCenterY(std::size_t i, std::size_t j, std::size_t k) const {
  return 0.25 * (_hz(i, j, k + 1) + _hz(i, j, k) + _hz(i, j - 1, k + 1) +
                 _hz(i, j - 1, k));
}

double EMF::hxFaceCenterY(std::size_t i, std::size_t j, std::size_t k) const {
  return 0.25 * (_hx(i + 1, j, k) + _hx(i, j, k) + _hx(i + 1, j - 1, k) +
                 _hx(i, j - 1, k));
}

double EMF::exFaceCenterZ(std::size_t i, std::size_t j, std::size_t k) const {
  return 0.5 * (_ex(i, j + 1, k) + _ex(i, j, k));
}

double EMF::eyFaceCenterZ(std::size_t i, std::size_t j, std::size_t k) const {
  return 0.5 * (_ey(i + 1, j, k) + _ey(i, j, k));
}

double EMF::hxFaceCenterZ(std::size_t i, std::size_t j, std::size_t k) const {
  return 0.25 * (_hx(i + 1, j, k) + _hx(i, j, k) + _hx(i + 1, j, k - 1) +
                 _hx(i, j, k - 1));
}

double EMF::hyFaceCenterZ(std::size_t i, std::size_t j, std::size_t k) const {
  return 0.25 * (_hy(i, j + 1, k) + _hy(i, j, k) + _hy(i, j + 1, k - 1) +
                 _hy(i, j, k - 1));
}

xt::xarray<double>& EMF::field(Attribute a, Component c) {
  return field(EMF::fieldFromAttributeAndComponent(a, c));
}

void EMF::allocateEx(std::size_t nx, std::size_t ny, std::size_t nz) {
  _ex = xt::zeros<double>({nx, ny, nz});
}

void EMF::allocateEy(std::size_t nx, std::size_t ny, std::size_t nz) {
  _ey = xt::zeros<double>({nx, ny, nz});
}

void EMF::allocateEz(std::size_t nx, std::size_t ny, std::size_t nz) {
  _ez = xt::zeros<double>({nx, ny, nz});
}

void EMF::allocateHx(std::size_t nx, std::size_t ny, std::size_t nz) {
  _hx = xt::zeros<double>({nx, ny, nz});
}

void EMF::allocateHy(std::size_t nx, std::size_t ny, std::size_t nz) {
  _hy = xt::zeros<double>({nx, ny, nz});
}

void EMF::allocateHz(std::size_t nx, std::size_t ny, std::size_t nz) {
  _hz = xt::zeros<double>({nx, ny, nz});
}

void EMF::allocateJx(std::size_t nx, std::size_t ny, std::size_t nz) {
  _jx = xt::zeros<double>({nx, ny, nz});
}

void EMF::allocateJy(std::size_t nx, std::size_t ny, std::size_t nz) {
  _jy = xt::zeros<double>({nx, ny, nz});
}

void EMF::allocateJz(std::size_t nx, std::size_t ny, std::size_t nz) {
  _jz = xt::zeros<double>({nx, ny, nz});
}

void EMF::allocateJxPrev(std::size_t nx, std::size_t ny, std::size_t nz) {
  _jx_prev = xt::zeros<double>({nx, ny, nz});
}

void EMF::allocateJyPrev(std::size_t nx, std::size_t ny, std::size_t nz) {
  _jy_prev = xt::zeros<double>({nx, ny, nz});
}

void EMF::allocateJzPrev(std::size_t nx, std::size_t ny, std::size_t nz) {
  _jz_prev = xt::zeros<double>({nx, ny, nz});
}

void EMF::allocateExPrev(std::size_t nx, std::size_t ny, std::size_t nz) {
  _ex_prev = xt::zeros<double>({nx, ny, nz});
}

void EMF::allocateEyPrev(std::size_t nx, std::size_t ny, std::size_t nz) {
  _ey_prev = xt::zeros<double>({nx, ny, nz});
}

void EMF::allocateEzPrev(std::size_t nx, std::size_t ny, std::size_t nz) {
  _ez_prev = xt::zeros<double>({nx, ny, nz});
}

}  // namespace xfdtd
