#ifndef _XFDTD_CORE_ELECTROMAGNETIC_FIELD_H_
#define _XFDTD_CORE_ELECTROMAGNETIC_FIELD_H_

#include <xfdtd/coordinate_system/coordinate_system.h>
#include <xfdtd/exception/exception.h>
#include <xfdtd/common/type_define.h>

#include <utility>

namespace xfdtd {

class XFDTDEMFException : public XFDTDException {
 public:
  explicit XFDTDEMFException(std::string message = "XFDTD EMF Exception")
      : XFDTDException(std::move(message)) {}
};

class EMF {
 public:
  enum class Component { X, Y, Z, Magnitude };

  enum class Attribute { E, H };

  enum class Field { EX, EY, EZ, EM, HX, HY, HZ, HM };

 public:
  template <EMF::Attribute a>
  static constexpr auto dualAttribute() -> Attribute;

 public:
  static Field fieldFromAttributeAndComponent(Attribute a, Component c);

  static Component componentFromField(Field f);

  static Attribute attributeFromField(Field f);

  const Array3D<Real>& ex() const;

  const Array3D<Real>& ey() const;

  const Array3D<Real>& ez() const;

  const Array3D<Real>& hx() const;

  const Array3D<Real>& hy() const;

  const Array3D<Real>& hz() const;

  const Array3D<Real>& field(Field f) const;

  const Array3D<Real>& field(Attribute a, Component c) const;

  Array3D<Real>& ex();

  Array3D<Real>& ey();

  Array3D<Real>& ez();

  Array3D<Real>& hx();

  Array3D<Real>& hy();

  Array3D<Real>& hz();

  Array3D<Real>& field(Field f);

  Array3D<Real>& field(Attribute a, Component c);

  void allocateEx(std::size_t nx, std::size_t ny, std::size_t nz);

  void allocateEy(std::size_t nx, std::size_t ny, std::size_t nz);

  void allocateEz(std::size_t nx, std::size_t ny, std::size_t nz);

  void allocateHx(std::size_t nx, std::size_t ny, std::size_t nz);

  void allocateHy(std::size_t nx, std::size_t ny, std::size_t nz);

  void allocateHz(std::size_t nx, std::size_t ny, std::size_t nz);

 private:
  Array3D<Real> _ex, _ey, _ez;
  Array3D<Real> _hx, _hy, _hz;
};

template <EMF::Attribute a>
inline constexpr auto EMF::dualAttribute() -> EMF::Attribute {
  if constexpr (a == EMF::Attribute::E) {
    return EMF::Attribute::H;
  } else if constexpr (a == EMF::Attribute::H) {
    return EMF::Attribute::E;
  } else {
    throw XFDTDEMFException{"dualAttribute(): Invalid EMF::Attribute"};
  }
}

template <EMF::Field f>
inline constexpr auto fieldToAttribute() -> EMF::Attribute {
  if constexpr (f == EMF::Field::EX || f == EMF::Field::EY ||
                f == EMF::Field::EZ) {
    return EMF::Attribute::E;
  } else if constexpr (f == EMF::Field::HX || f == EMF::Field::HY ||
                       f == EMF::Field::HZ) {
    return EMF::Attribute::H;
  } else {
    throw XFDTDEMFException{"fieldToAttribute(): Invalid EMF::Field"};
  }
}

template <EMF::Field f>
inline constexpr auto fieldToComponent() -> EMF::Component {
  if constexpr (f == EMF::Field::EX || f == EMF::Field::HX) {
    return EMF::Component::X;
  } else if constexpr (f == EMF::Field::EY || f == EMF::Field::HY) {
    return EMF::Component::Y;
  } else if constexpr (f == EMF::Field::EZ || f == EMF::Field::HZ) {
    return EMF::Component::Z;
  } else if constexpr (f == EMF::Field::EM || f == EMF::Field::HM) {
    return EMF::Component::Magnitude;
  } else {
    throw XFDTDEMFException{"fieldToComponent(): Invalid EMF::Field"};
  }
}

template <EMF::Attribute a, EMF::Component c>
inline constexpr auto attributeComponentToField() -> EMF::Field {
  if constexpr (a == EMF::Attribute::E) {
    if constexpr (c == EMF::Component::X) {
      return EMF::Field::EX;
    } else if constexpr (c == EMF::Component::Y) {
      return EMF::Field::EY;
    } else if constexpr (c == EMF::Component::Z) {
      return EMF::Field::EZ;
    } else if constexpr (c == EMF::Component::Magnitude) {
      return EMF::Field::EM;
    } else {
      throw XFDTDEMFException{
          "attributeComponentToField(): Invalid EMF::Component"};
    }
  } else if constexpr (a == EMF::Attribute::H) {
    if constexpr (c == EMF::Component::X) {
      return EMF::Field::HX;
    } else if constexpr (c == EMF::Component::Y) {
      return EMF::Field::HY;
    } else if constexpr (c == EMF::Component::Z) {
      return EMF::Field::HZ;
    } else if constexpr (c == EMF::Component::Magnitude) {
      return EMF::Field::HM;
    } else {
      throw XFDTDEMFException{
          "attributeComponentToField(): Invalid EMF::Component"};
    }
  } else {
    throw XFDTDEMFException{
        "attributeComponentToField(): Invalid EMF::Attribute"};
  }
}

}  // namespace xfdtd

#endif  // _XFDTD_CORE_ELECTROMAGNETIC_FIELD_H_
