#include <xfdtd/electromagnetic_field/electromagnetic_field.h>

#include <cstddef>

namespace xfdtd {

const Array3D<Real>& EMF::ex() const { return _ex; }

const Array3D<Real>& EMF::ey() const { return _ey; }

const Array3D<Real>& EMF::ez() const { return _ez; }

const Array3D<Real>& EMF::hx() const { return _hx; }

const Array3D<Real>& EMF::hy() const { return _hy; }

const Array3D<Real>& EMF::hz() const { return _hz; }

const Array3D<Real>& EMF::field(Field f) const {
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
    default:
      throw XFDTDEMFException("Invalid field type");
  }
}

const Array3D<Real>& EMF::field(Attribute a, Component c) const {
  return field(EMF::attributeComponentToField(a, c));
}

Array3D<Real>& EMF::ex() { return _ex; }

Array3D<Real>& EMF::ey() { return _ey; }

Array3D<Real>& EMF::ez() { return _ez; }

Array3D<Real>& EMF::hx() { return _hx; }

Array3D<Real>& EMF::hy() { return _hy; }

Array3D<Real>& EMF::hz() { return _hz; }

Array3D<Real>& EMF::field(Field f) {
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
    default:
      throw XFDTDEMFException("Invalid field type");
  }
}

Array3D<Real>& EMF::field(Attribute a, Component c) {
  return field(EMF::attributeComponentToField(a, c));
}

void EMF::allocateEx(std::size_t nx, std::size_t ny, std::size_t nz) {
  _ex = xt::zeros<Real>({nx, ny, nz});
}

void EMF::allocateEy(std::size_t nx, std::size_t ny, std::size_t nz) {
  _ey = xt::zeros<Real>({nx, ny, nz});
}

void EMF::allocateEz(std::size_t nx, std::size_t ny, std::size_t nz) {
  _ez = xt::zeros<Real>({nx, ny, nz});
}

void EMF::allocateHx(std::size_t nx, std::size_t ny, std::size_t nz) {
  _hx = xt::zeros<Real>({nx, ny, nz});
}

void EMF::allocateHy(std::size_t nx, std::size_t ny, std::size_t nz) {
  _hy = xt::zeros<Real>({nx, ny, nz});
}

void EMF::allocateHz(std::size_t nx, std::size_t ny, std::size_t nz) {
  _hz = xt::zeros<Real>({nx, ny, nz});
}

}  // namespace xfdtd
