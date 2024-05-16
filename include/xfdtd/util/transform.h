#ifndef _XFDTD_CORE_TRANSFORM_H_
#define _XFDTD_CORE_TRANSFORM_H_

#include <xfdtd/coordinate_system/coordinate_system.h>
#include <xfdtd/electromagnetic_field/electromagnetic_field.h>
#include <xfdtd/exception/exception.h>

#include <tuple>

namespace xfdtd::transform {

class XFDTDTransformException : public XFDTDException {
 public:
  explicit XFDTDTransformException(
      std::string message = "XFDTD Transform Exception")
      : XFDTDException(std::move(message)) {}
};

enum class SCS {
  R,
  Theta,
  Phi,
};

template <typename T>
inline auto xYZToABC(const std::tuple<T, T, T> &data, Axis::XYZ xyz) {
  switch (xyz) {
    case Axis::XYZ::X: {
      // a: y b: z c: x
      return std::make_tuple(std::get<1>(data), std::get<2>(data),
                             std::get<0>(data));
    }
    case Axis::XYZ::Y: {
      // a: z b: x c: y
      return std::make_tuple(std::get<2>(data), std::get<0>(data),
                             std::get<1>(data));
    }
    case Axis::XYZ::Z: {
      // a: x b: y c: z
      return std::make_tuple(std::get<0>(data), std::get<1>(data),
                             std::get<2>(data));
    }
    default:
      throw std::invalid_argument("Invalid Axis::XYZ");
  }
}

template <typename T>
inline auto aBCToXYZ(const std::tuple<T, T, T> &data, Axis::XYZ xyz) {
  switch (xyz) {
    case Axis::XYZ::X: {
      // x: c y: a z: b
      return std::make_tuple(std::get<2>(data), std::get<0>(data),
                             std::get<1>(data));
    }
    case Axis::XYZ::Y: {
      // x: b y: c z: a
      return std::make_tuple(std::get<1>(data), std::get<2>(data),
                             std::get<0>(data));
    }
    case Axis::XYZ::Z: {
      // x: a y: b z: c
      return std::make_tuple(std::get<0>(data), std::get<1>(data),
                             std::get<2>(data));
    }
    default:
      throw std::invalid_argument("Invalid Axis::XYZ");
  }
}

template <Axis::XYZ xyz, SCS scs>
inline auto cCSToSCSTransformMatrix(const auto &theta, const auto &phi) {
  if constexpr (scs == SCS::R) {
    if constexpr (xyz == Axis::XYZ::X) {
      return xt::sin(theta) * xt::cos(phi);
    } else if constexpr (xyz == Axis::XYZ::Y) {
      return xt::sin(theta) * xt::sin(phi);
    } else if constexpr (xyz == Axis::XYZ::Z) {
      return xt::cos(theta);
    } else {
      throw XFDTDException{"Invalid Axis::XYZ"};
    }
  } else if constexpr (scs == SCS::Theta) {
    if constexpr (xyz == Axis::XYZ::X) {
      return xt::cos(theta) * xt::cos(phi);
    } else if constexpr (xyz == Axis::XYZ::Y) {
      return xt::cos(theta) * xt::sin(phi);
    } else if constexpr (xyz == Axis::XYZ::Z) {
      return -1 * xt::sin(theta);
    } else {
      throw XFDTDException{"Invalid Axis::XYZ"};
    }
  } else if constexpr (scs == SCS::Phi) {
    if constexpr (xyz == Axis::XYZ::X) {
      return -1 * xt::sin(phi);
    } else if constexpr (xyz == Axis::XYZ::Y) {
      return xt::cos(phi);
    } else if constexpr (xyz == Axis::XYZ::Z) {
      return xt::zeros_like(theta);
    } else {
      throw XFDTDException{"Invalid Axis::XYZ"};
    }
  } else {
    throw XFDTDException{"Invalid SCS"};
  }
}

// Warning: This function is not tested.
template <Axis::XYZ xyz, SCS scs>
inline auto cCSToSCS(const auto &theta, const auto &phi, const auto &data) {
  if constexpr (scs == SCS::R) {
    if constexpr (xyz == Axis::XYZ::X) {
      auto &&sin_t_cos_p = xt::sin(theta) * xt::cos(phi);
      return sin_t_cos_p * data;
    } else if constexpr (xyz == Axis::XYZ::Y) {
      auto &&sin_t_sin_p = xt::sin(theta) * xt::sin(phi);
      return sin_t_sin_p * data;
    } else if constexpr (xyz == Axis::XYZ::Z) {
      auto &&cos_t = xt::cos(theta);
      return cos_t * data;
    } else {
      throw XFDTDException{"Invalid Axis::XYZ"};
    }
  } else if constexpr (scs == SCS::Theta) {
    if constexpr (xyz == Axis::XYZ::X) {
      auto &&cos_t_cos_p = xt::cos(theta) * xt::cos(phi);
      return cos_t_cos_p * data;
    } else if constexpr (xyz == Axis::XYZ::Y) {
      auto &&cos_t_sin_p = xt::cos(theta) * xt::sin(phi);
      return cos_t_sin_p * data;
    } else if constexpr (xyz == Axis::XYZ::Z) {
      auto &&n_sin_t = -1 * xt::sin(theta);
      return n_sin_t * data;
    } else {
      throw XFDTDException{"Invalid Axis::XYZ"};
    }
  } else if constexpr (scs == SCS::Phi) {
    if constexpr (xyz == Axis::XYZ::X) {
      auto &&n_sin_p = -1 * xt::sin(phi);
      return n_sin_p * data;
    } else if constexpr (xyz == Axis::XYZ::Y) {
      auto &&cos_p = xt::cos(phi);
      return cos_p * data;
    } else if constexpr (xyz == Axis::XYZ::Z) {
      return xt::empty_like(data);
    } else {
      throw XFDTDException{"Invalid Axis::XYZ"};
    }
  } else {
    throw XFDTDException{"Invalid SCS"};
  }
}

template <Axis::XYZ xyz>
inline constexpr auto xYZToComponent() -> EMF::Component {
  if constexpr (xyz == Axis::XYZ::X) {
    return EMF::Component::X;
  } else if constexpr (xyz == Axis::XYZ::Y) {
    return EMF::Component::Y;
  } else if constexpr (xyz == Axis::XYZ::Z) {
    return EMF::Component::Z;
  } else {
    throw XFDTDTransformException{"fromXYZToComponent(): Invalid Axis::XYZ"};
  }
}

template <EMF::Attribute a, Axis::XYZ xyz>
inline constexpr auto attributeXYZToField() -> EMF::Field {
  if constexpr (a == EMF::Attribute::E) {
    if constexpr (xyz == Axis::XYZ::X) {
      return EMF::Field::EX;
    } else if constexpr (xyz == Axis::XYZ::Y) {
      return EMF::Field::EY;
    } else if constexpr (xyz == Axis::XYZ::Z) {
      return EMF::Field::EZ;
    } else {
      throw XFDTDTransformException{"attributeXYZToField(): Invalid Axis::XYZ"};
    }
  } else if constexpr (a == EMF::Attribute::H) {
    if constexpr (xyz == Axis::XYZ::X) {
      return EMF::Field::HX;
    } else if constexpr (xyz == Axis::XYZ::Y) {
      return EMF::Field::HY;
    } else if constexpr (xyz == Axis::XYZ::Z) {
      return EMF::Field::HZ;
    } else {
      throw XFDTDTransformException{"attributeXYZToField(): Invalid Axis::XYZ"};
    }
  } else {
    throw XFDTDTransformException{
        "attributeXYZToField(): Invalid EMF::Attribute"};
  }
}

template <EMF::Attribute a, Axis::XYZ xyz>
inline constexpr auto attributeXYZToTangentialAAxis() -> EMF::Field {
  if constexpr (a == EMF::Attribute::E) {
    return attributeXYZToField<a, Axis::tangentialAAxis<xyz>>();
  } else {
    throw XFDTDTransformException{
        "attributeXYZToTangentialAAxis(): Invalid EMF::Attribute"};
  }
}

}  // namespace xfdtd::transform

#endif  // _XFDTD_CORE_TRANSFORM_H_
