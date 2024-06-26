#ifndef _XFDTD_CORE_ELECTROMAGNETIC_FIELD_H_
#define _XFDTD_CORE_ELECTROMAGNETIC_FIELD_H_

#include <xfdtd/common/type_define.h>
#include <xfdtd/coordinate_system/coordinate_system.h>
#include <xfdtd/exception/exception.h>

#include <utility>

namespace xfdtd {

class XFDTDEMFException : public XFDTDException {
 public:
  explicit XFDTDEMFException(std::string message)
      : XFDTDException(std::move(message)) {}
};

class EMF {
 public:
  enum class Component { X, Y, Z, Magnitude };

  enum class Attribute { E, H };

  enum class Field { EX, EY, EZ, EM, HX, HY, HZ, HM };

 public:
  static constexpr auto componentToString(Component c) -> std::string;

  static constexpr auto attributeToString(Attribute a) -> std::string;

  static constexpr auto fieldToString(Field f) -> std::string;

  static constexpr auto attributeComponentToField(Attribute a,
                                                  Component c) -> Field;

  static constexpr auto fieldToComponent(Field f) -> Component;

  static constexpr auto fieldToAttribute(Field f) -> Attribute;

  static constexpr auto dualAttribute(Attribute a) -> Attribute;

  static constexpr auto xYZToComponent(Axis::XYZ xyz) -> Component;

  static constexpr auto componentToXYZ(Component c) -> Axis::XYZ;

 public:
  template <Field f>
  auto field() const -> const Array3D<Real>&;

  template <Field f>
  auto field() -> Array3D<Real>&;

 public:
  template <Attribute a, Axis::XYZ xyz>
  auto field() const -> const Array3D<Real>& {
    constexpr auto f = attributeComponentToField(a, xYZToComponent(xyz));
    if constexpr (f == Field::EX) {
      return _ex;
    } else if constexpr (f == Field::EY) {
      return _ey;
    } else if constexpr (f == Field::EZ) {
      return _ez;
    } else if constexpr (f == Field::HX) {
      return _hx;
    } else if constexpr (f == Field::HY) {
      return _hy;
    } else if constexpr (f == Field::HZ) {
      return _hz;
    }
  }

  template <Attribute a, Axis::XYZ xyz>
  auto field() -> Array3D<Real>& {
    return const_cast<Array3D<Real>&>(
        static_cast<const EMF*>(this)->field<a, xyz>());
  }

 public:
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

// inline constexpr auto EMF::componentToString(Component c) -> std::string {
//   // is it possible to use switch case here?
//   if (c == Component::X) {
//     return "X";
//   }
//   if (c == Component::Y) {
//     return "Y";
//   }
//   if (c == Component::Z) {
//     return "Z";
//   }
//   if (c == Component::Magnitude) {
//     return "Magnitude";
//   }

//   return "Invalid Component";
// }

// inline constexpr auto EMF::attributeToString(Attribute a) -> std::string {
//   if (a == Attribute::E) {
//     return "E";
//   }
//   if (a == Attribute::H) {
//     return "H";
//   }

//   return "Invalid Attribute";
// }

// inline constexpr auto EMF::fieldToString(Field f) -> std::string {
//   if (f == Field::EX) {
//     return "Ex";
//   }
//   if (f == Field::EY) {
//     return "Ey";
//   }
//   if (f == Field::EZ) {
//     return "Ez";
//   }
//   if (f == Field::EM) {
//     return "Em";
//   }
//   if (f == Field::HX) {
//     return "Hx";
//   }
//   if (f == Field::HY) {
//     return "Hy";
//   }
//   if (f == Field::HZ) {
//     return "Hz";
//   }
//   if (f == Field::HM) {
//     return "Hm";
//   }

//   return "Invalid Field";
// }

inline constexpr auto EMF::attributeComponentToField(
    EMF::Attribute a, EMF::Component c) -> EMF::Field {
  if (a == EMF::Attribute::E) {
    if (c == EMF::Component::X) {
      return EMF::Field::EX;
    }
    if (c == EMF::Component::Y) {
      return EMF::Field::EY;
    }
    if (c == EMF::Component::Z) {
      return EMF::Field::EZ;
    }
    if (c == EMF::Component::Magnitude) {
      return EMF::Field::EM;
    }
  }

  if (a == EMF::Attribute::H) {
    if (c == EMF::Component::X) {
      return EMF::Field::HX;
    }
    if (c == EMF::Component::Y) {
      return EMF::Field::HY;
    }
    if (c == EMF::Component::Z) {
      return EMF::Field::HZ;
    }
    if (c == EMF::Component::Magnitude) {
      return EMF::Field::HM;
    }
  }

  throw XFDTDEMFException("Invalid Attribute or Component");
}

inline constexpr EMF::Component EMF::fieldToComponent(EMF::Field f) {
  if (f == EMF::Field::EX || f == EMF::Field::HX) {
    return EMF::Component::X;
  }
  if (f == EMF::Field::EY || f == EMF::Field::HY) {
    return EMF::Component::Y;
  }
  if (f == EMF::Field::EZ || f == EMF::Field::HZ) {
    return EMF::Component::Z;
  }
  if (f == EMF::Field::EM || f == EMF::Field::HM) {
    return EMF::Component::Magnitude;
  }

  throw XFDTDEMFException("Invalid Field");
}

constexpr EMF::Attribute EMF::fieldToAttribute(EMF::Field f) {
  if (f == EMF::Field::EX || f == EMF::Field::EY || f == EMF::Field::EZ ||
      f == EMF::Field::EM) {
    return EMF::Attribute::E;
  }
  if (f == EMF::Field::HX || f == EMF::Field::HY || f == EMF::Field::HZ ||
      f == EMF::Field::HM) {
    return EMF::Attribute::H;
  }

  throw XFDTDEMFException("Invalid Field");
}

inline constexpr auto EMF::dualAttribute(EMF::Attribute a) -> EMF::Attribute {
  if (a == EMF::Attribute::E) {
    return EMF::Attribute::H;
  }
  if (a == EMF::Attribute::H) {
    return EMF::Attribute::E;
  }

  throw XFDTDEMFException("Invalid Attribute");
}

inline constexpr auto EMF::xYZToComponent(Axis::XYZ xyz) -> EMF::Component {
  if (xyz == Axis::XYZ::X) {
    return EMF::Component::X;
  }
  if (xyz == Axis::XYZ::Y) {
    return EMF::Component::Y;
  }
  if (xyz == Axis::XYZ::Z) {
    return EMF::Component::Z;
  }

  throw XFDTDEMFException("xYZToComponent(): Invalid Axis::XYZ");
}

inline constexpr auto EMF::componentToXYZ(EMF::Component c) -> Axis::XYZ {
  if (c == EMF::Component::X) {
    return Axis::XYZ::X;
  }
  if (c == EMF::Component::Y) {
    return Axis::XYZ::Y;
  }
  if (c == EMF::Component::Z) {
    return Axis::XYZ::Z;
  }

  throw XFDTDEMFException("componentToXYZ(): Invalid EMF::Component");
}

template <EMF::Field f>
inline auto EMF::field() const -> const Array3D<Real>& {
  // static assertion to check if the function is called with valid template
  // parameter
  static_assert(
      []() {
        return f == EMF::Field::EX || f == EMF::Field::EY ||
               f == EMF::Field::EZ || f == EMF::Field::HX ||
               f == EMF::Field::HY || f == EMF::Field::HZ;
      }(),
      "EMF::field(): Invalid EMF::Field");

  if constexpr (f == EMF::Field::EX) {
    return _ex;
  } else if constexpr (f == EMF::Field::EY) {
    return _ey;
  } else if constexpr (f == EMF::Field::EZ) {
    return _ez;
  } else if constexpr (f == EMF::Field::HX) {
    return _hx;
  } else if constexpr (f == EMF::Field::HY) {
    return _hy;
  } else if constexpr (f == EMF::Field::HZ) {
    return _hz;
  } else {
    throw XFDTDEMFException{"field(): Invalid EMF::Field"};
  }
}

template <EMF::Field f>
inline auto EMF::field() -> Array3D<Real>& {
  return const_cast<Array3D<Real>&>(static_cast<const EMF*>(this)->field<f>());
}

}  // namespace xfdtd

#endif  // _XFDTD_CORE_ELECTROMAGNETIC_FIELD_H_
