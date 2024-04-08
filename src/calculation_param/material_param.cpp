#include <xfdtd/calculation_param/calculation_param.h>
#include <xfdtd/coordinate_system/coordinate_system.h>
#include <xfdtd/electromagnetic_field/electromagnetic_field.h>
#include <xfdtd/material/material.h>
#include <xfdtd/common/constant.h>
#include <xfdtd/util/fdtd_basic.h>

#include <memory>
#include <utility>
#include <vector>

namespace xfdtd {

MaterialParam::Attribute MaterialParam::fromPropertyToAttribute(
    Property property) {
  switch (property) {
    case MaterialParam::Property::EPS_X:
    case MaterialParam::Property::EPS_Y:
    case MaterialParam::Property::EPS_Z:
      return MaterialParam::Attribute::EPSILON;
    case MaterialParam::Property::MU_X:
    case MaterialParam::Property::MU_Y:
    case MaterialParam::Property::MU_Z:
      return MaterialParam::Attribute::MU;
    case MaterialParam::Property::SIGMA_E_X:
    case MaterialParam::Property::SIGMA_E_Y:
    case MaterialParam::Property::SIGMA_E_Z:
      return MaterialParam::Attribute::SIGMA_E;
    case MaterialParam::Property::SIGMA_M_X:
    case MaterialParam::Property::SIGMA_M_Y:
    case MaterialParam::Property::SIGMA_M_Z:
      return MaterialParam::Attribute::SIGMA_M;
    default:
      throw XFDTDCalculationParamException("Invalid property");
  }
}

Axis::XYZ MaterialParam::fromPropertyToXYZ(Property property) {
  switch (property) {
    case MaterialParam::Property::EPS_X:
    case MaterialParam::Property::MU_X:
    case MaterialParam::Property::SIGMA_E_X:
    case MaterialParam::Property::SIGMA_M_X:
      return Axis::XYZ::X;
    case MaterialParam::Property::EPS_Y:
    case MaterialParam::Property::MU_Y:
    case MaterialParam::Property::SIGMA_E_Y:
    case MaterialParam::Property::SIGMA_M_Y:
      return Axis::XYZ::Y;
    case MaterialParam::Property::EPS_Z:
    case MaterialParam::Property::MU_Z:
    case MaterialParam::Property::SIGMA_E_Z:
    case MaterialParam::Property::SIGMA_M_Z:
      return Axis::XYZ::Z;
    default:
      throw XFDTDCalculationParamException("Invalid property");
  }
}

MaterialParam::Property MaterialParam::fromAttributeAndXYZ(Attribute attribute,
                                                           Axis::XYZ xyz) {
  switch (attribute) {
    case Attribute::EPSILON:
      switch (xyz) {
        case Axis::XYZ::X:
          return Property::EPS_X;
        case Axis::XYZ::Y:
          return Property::EPS_Y;
        case Axis::XYZ::Z:
          return Property::EPS_Z;
        default:
          throw XFDTDCalculationParamException("Invalid xyz");
      }
    case Attribute::MU:
      switch (xyz) {
        case Axis::XYZ::X:
          return Property::MU_X;
        case Axis::XYZ::Y:
          return Property::MU_Y;
        case Axis::XYZ::Z:
          return Property::MU_Z;
        default:
          throw XFDTDCalculationParamException("Invalid xyz");
      }
    case Attribute::SIGMA_E:
      switch (xyz) {
        case Axis::XYZ::X:
          return Property::SIGMA_E_X;
        case Axis::XYZ::Y:
          return Property::SIGMA_E_Y;
        case Axis::XYZ::Z:
          return Property::SIGMA_E_Z;
        default:
          throw XFDTDCalculationParamException("Invalid xyz");
      }
    case Attribute::SIGMA_M:
      switch (xyz) {
        case Axis::XYZ::X:
          return Property::SIGMA_M_X;
        case Axis::XYZ::Y:
          return Property::SIGMA_M_Y;
        case Axis::XYZ::Z:
          return Property::SIGMA_M_Z;
        default:
          throw XFDTDCalculationParamException("Invalid xyz");
      }
    default:
      throw XFDTDCalculationParamException("Invalid attribute");
  }
}

const Array3D<Real>& MaterialParam::epsX() const { return _eps_x; }

const Array3D<Real>& MaterialParam::epsY() const { return _eps_y; }

const Array3D<Real>& MaterialParam::epsZ() const { return _eps_z; }

const Array3D<Real>& MaterialParam::muX() const { return _mu_x; }

const Array3D<Real>& MaterialParam::muY() const { return _mu_y; }

const Array3D<Real>& MaterialParam::muZ() const { return _mu_z; }

const Array3D<Real>& MaterialParam::sigmaEX() const { return _sigma_e_x; }

const Array3D<Real>& MaterialParam::sigmaEY() const { return _sigma_e_y; }

const Array3D<Real>& MaterialParam::sigmaEZ() const { return _sigma_e_z; }

const Array3D<Real>& MaterialParam::sigmaMX() const { return _sigma_m_x; }

const Array3D<Real>& MaterialParam::sigmaMY() const { return _sigma_m_y; }

const Array3D<Real>& MaterialParam::sigmaMZ() const { return _sigma_m_z; }

const Array3D<Real>& MaterialParam::property(
    MaterialParam::Property property) const {
  switch (property) {
    case MaterialParam::Property::EPS_X:
      return _eps_x;
    case MaterialParam::Property::EPS_Y:
      return _eps_y;
    case MaterialParam::Property::EPS_Z:
      return _eps_z;
    case MaterialParam::Property::MU_X:
      return _mu_x;
    case MaterialParam::Property::MU_Y:
      return _mu_y;
    case MaterialParam::Property::MU_Z:
      return _mu_z;
    case MaterialParam::Property::SIGMA_E_X:
      return _sigma_e_x;
    case MaterialParam::Property::SIGMA_E_Y:
      return _sigma_e_y;
    case MaterialParam::Property::SIGMA_E_Z:
      return _sigma_e_z;
    case MaterialParam::Property::SIGMA_M_X:
      return _sigma_m_x;
    case MaterialParam::Property::SIGMA_M_Y:
      return _sigma_m_y;
    case MaterialParam::Property::SIGMA_M_Z:
      return _sigma_m_z;
    default:
      throw XFDTDCalculationParamException("Invalid property");
  }
}

const Array3D<Real>& MaterialParam::property(MaterialParam::Attribute attribute,
                                             Axis::XYZ xyz) const {
  switch (attribute) {
    case MaterialParam::Attribute::EPSILON:
      switch (xyz) {
        case Axis::XYZ::X:
          return _eps_x;
        case Axis::XYZ::Y:
          return _eps_y;
        case Axis::XYZ::Z:
          return _eps_z;
        default:
          throw XFDTDCalculationParamException("Invalid xyz");
      }
    case MaterialParam::Attribute::MU:
      switch (xyz) {
        case Axis::XYZ::X:
          return _mu_x;
        case Axis::XYZ::Y:
          return _mu_y;
        case Axis::XYZ::Z:
          return _mu_z;
        default:
          throw XFDTDCalculationParamException("Invalid xyz");
      }
    case MaterialParam::Attribute::SIGMA_E:
      switch (xyz) {
        case Axis::XYZ::X:
          return _sigma_e_x;
        case Axis::XYZ::Y:
          return _sigma_e_y;
        case Axis::XYZ::Z:
          return _sigma_e_z;
        default:
          throw XFDTDCalculationParamException("Invalid xyz");
      }
    case MaterialParam::Attribute::SIGMA_M:
      switch (xyz) {
        case Axis::XYZ::X:
          return _sigma_m_x;
        case Axis::XYZ::Y:
          return _sigma_m_y;
        case Axis::XYZ::Z:
          return _sigma_m_z;
        default:
          throw XFDTDCalculationParamException("Invalid xyz");
      }
    default:
      throw XFDTDCalculationParamException("Invalid attribute");
  }
}

void MaterialParam::addMaterial(std::shared_ptr<Material> material) {
  _materials.emplace_back(std::move(material));
}

Array3D<Real>& MaterialParam::epsX() { return _eps_x; }

Array3D<Real>& MaterialParam::epsY() { return _eps_y; }

Array3D<Real>& MaterialParam::epsZ() { return _eps_z; }

Array3D<Real>& MaterialParam::muX() { return _mu_x; }

Array3D<Real>& MaterialParam::muY() { return _mu_y; }

Array3D<Real>& MaterialParam::muZ() { return _mu_z; }

Array3D<Real>& MaterialParam::sigmaEX() { return _sigma_e_x; }

Array3D<Real>& MaterialParam::sigmaEY() { return _sigma_e_y; }

Array3D<Real>& MaterialParam::sigmaEZ() { return _sigma_e_z; }

Array3D<Real>& MaterialParam::sigmaMX() { return _sigma_m_x; }

Array3D<Real>& MaterialParam::sigmaMY() { return _sigma_m_y; }

Array3D<Real>& MaterialParam::sigmaMZ() { return _sigma_m_z; }

Array3D<Real>& MaterialParam::property(MaterialParam::Property property) {
  switch (property) {
    case MaterialParam::Property::EPS_X:
      return _eps_x;
    case MaterialParam::Property::EPS_Y:
      return _eps_y;
    case MaterialParam::Property::EPS_Z:
      return _eps_z;
    case MaterialParam::Property::MU_X:
      return _mu_x;
    case MaterialParam::Property::MU_Y:
      return _mu_y;
    case MaterialParam::Property::MU_Z:
      return _mu_z;
    case MaterialParam::Property::SIGMA_E_X:
      return _sigma_e_x;
    case MaterialParam::Property::SIGMA_E_Y:
      return _sigma_e_y;
    case MaterialParam::Property::SIGMA_E_Z:
      return _sigma_e_z;
    case MaterialParam::Property::SIGMA_M_X:
      return _sigma_m_x;
    case MaterialParam::Property::SIGMA_M_Y:
      return _sigma_m_y;
    case MaterialParam::Property::SIGMA_M_Z:
      return _sigma_m_z;
    default:
      throw XFDTDCalculationParamException("Invalid property");
  }
}

Array3D<Real>& MaterialParam::property(MaterialParam::Attribute attribute,
                                       Axis::XYZ xyz) {
  switch (attribute) {
    case MaterialParam::Attribute::EPSILON:
      switch (xyz) {
        case Axis::XYZ::X:
          return _eps_x;
        case Axis::XYZ::Y:
          return _eps_y;
        case Axis::XYZ::Z:
          return _eps_z;
        default:
          throw XFDTDCalculationParamException("Invalid xyz");
      }
    case MaterialParam::Attribute::MU:
      switch (xyz) {
        case Axis::XYZ::X:
          return _mu_x;
        case Axis::XYZ::Y:
          return _mu_y;
        case Axis::XYZ::Z:
          return _mu_z;
        default:
          throw XFDTDCalculationParamException("Invalid xyz");
      }
    case MaterialParam::Attribute::SIGMA_E:
      switch (xyz) {
        case Axis::XYZ::X:
          return _sigma_e_x;
        case Axis::XYZ::Y:
          return _sigma_e_y;
        case Axis::XYZ::Z:
          return _sigma_e_z;
        default:
          throw XFDTDCalculationParamException("Invalid xyz");
      }
    case MaterialParam::Attribute::SIGMA_M:
      switch (xyz) {
        case Axis::XYZ::X:
          return _sigma_m_x;
        case Axis::XYZ::Y:
          return _sigma_m_y;
        case Axis::XYZ::Z:
          return _sigma_m_z;
        default:
          throw XFDTDCalculationParamException("Invalid xyz");
      }
    default:
      throw XFDTDCalculationParamException("Invalid attribute");
  }
}

template <EMF::Attribute a, Axis::XYZ xyz>
inline static auto makeArr(std::size_t nx, std::size_t ny, std::size_t nz)
    -> Array3D<Real> {
  if constexpr (a == EMF::Attribute::E) {
    if constexpr (xyz == Axis::XYZ::X) {
      const auto e_x_n = basic::GridStructure::exSize(nx, ny, nz);
      return xt::zeros<Real>(e_x_n);
    } else if constexpr (xyz == Axis::XYZ::Y) {
      const auto e_y_n = basic::GridStructure::eySize(nx, ny, nz);
      return xt::zeros<Real>(e_y_n);
    } else if constexpr (xyz == Axis::XYZ::Z) {
      const auto e_z_n = basic::GridStructure::ezSize(nx, ny, nz);
      return xt::zeros<Real>(e_z_n);
    } else {
      throw XFDTDCalculationParamException{"Invalid Axis::XYZ"};
    }
  } else if constexpr (a == EMF::Attribute::H) {
    if constexpr (xyz == Axis::XYZ::X) {
      const auto h_x_n = basic::GridStructure::hxSize(nx, ny, nz);
      return xt::zeros<Real>(h_x_n);
    } else if constexpr (xyz == Axis::XYZ::Y) {
      const auto h_y_n = basic::GridStructure::hySize(nx, ny, nz);
      return xt::zeros<Real>(h_y_n);
    } else if constexpr (xyz == Axis::XYZ::Z) {
      const auto h_z_n = basic::GridStructure::hzSize(nx, ny, nz);
      return xt::zeros<Real>(h_z_n);
    } else {
      throw XFDTDCalculationParamException{"Invalid Axis::XYZ"};
    }
  } else {
    throw XFDTDCalculationParamException{"Invalid EMF::Attribute"};
  }
}

void MaterialParam::allocate(std::size_t nx, std::size_t ny, std::size_t nz) {
  _eps_x = makeArr<EMF::Attribute::E, Axis::XYZ::X>(nx, ny, nz);
  _eps_x.fill(constant::EPSILON_0);
  _eps_y = makeArr<EMF::Attribute::E, Axis::XYZ::Y>(nx, ny, nz);
  _eps_y.fill(constant::EPSILON_0);
  _eps_z = makeArr<EMF::Attribute::E, Axis::XYZ::Z>(nx, ny, nz);
  _eps_z.fill(constant::EPSILON_0);

  _mu_x = makeArr<EMF::Attribute::H, Axis::XYZ::X>(nx, ny, nz);
  _mu_x.fill(constant::MU_0);
  _mu_y = makeArr<EMF::Attribute::H, Axis::XYZ::Y>(nx, ny, nz);
  _mu_y.fill(constant::MU_0);
  _mu_z = makeArr<EMF::Attribute::H, Axis::XYZ::Z>(nx, ny, nz);
  _mu_z.fill(constant::MU_0);

  _sigma_e_x = makeArr<EMF::Attribute::E, Axis::XYZ::X>(nx, ny, nz);
  _sigma_e_x.fill(constant::SIGMA_E_ZERO_APPROX);
  _sigma_e_y = makeArr<EMF::Attribute::E, Axis::XYZ::Y>(nx, ny, nz);
  _sigma_e_y.fill(constant::SIGMA_E_ZERO_APPROX);
  _sigma_e_z = makeArr<EMF::Attribute::E, Axis::XYZ::Z>(nx, ny, nz);
  _sigma_e_z.fill(constant::SIGMA_E_ZERO_APPROX);

  _sigma_m_x = makeArr<EMF::Attribute::H, Axis::XYZ::X>(nx, ny, nz);
  _sigma_m_x.fill(constant::SIGMA_M_ZERO_APPROX);
  _sigma_m_y = makeArr<EMF::Attribute::H, Axis::XYZ::Y>(nx, ny, nz);
  _sigma_m_y.fill(constant::SIGMA_M_ZERO_APPROX);
  _sigma_m_z = makeArr<EMF::Attribute::H, Axis::XYZ::Z>(nx, ny, nz);
  _sigma_m_z.fill(constant::SIGMA_M_ZERO_APPROX);
}

}  // namespace xfdtd
