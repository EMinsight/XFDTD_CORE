#include <xfdtd/calculation_param/calculation_param.h>

#include <memory>
#include <utility>
#include <vector>

#include "xfdtd/coordinate_system/coordinate_system.h"
#include "xfdtd/material/material.h"

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

const xt::xarray<double>& MaterialParam::epsX() const { return _eps_x; }

const xt::xarray<double>& MaterialParam::epsY() const { return _eps_y; }

const xt::xarray<double>& MaterialParam::epsZ() const { return _eps_z; }

const xt::xarray<double>& MaterialParam::muX() const { return _mu_x; }

const xt::xarray<double>& MaterialParam::muY() const { return _mu_y; }

const xt::xarray<double>& MaterialParam::muZ() const { return _mu_z; }

const xt::xarray<double>& MaterialParam::sigmaEX() const { return _sigma_e_x; }

const xt::xarray<double>& MaterialParam::sigmaEY() const { return _sigma_e_y; }

const xt::xarray<double>& MaterialParam::sigmaEZ() const { return _sigma_e_z; }

const xt::xarray<double>& MaterialParam::sigmaMX() const { return _sigma_m_x; }

const xt::xarray<double>& MaterialParam::sigmaMY() const { return _sigma_m_y; }

const xt::xarray<double>& MaterialParam::sigmaMZ() const { return _sigma_m_z; }

const xt::xarray<double>& MaterialParam::property(
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

const xt::xarray<double>& MaterialParam::property(
    MaterialParam::Attribute attribute, Axis::XYZ xyz) const {
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

xt::xarray<double>& MaterialParam::epsX() { return _eps_x; }

xt::xarray<double>& MaterialParam::epsY() { return _eps_y; }

xt::xarray<double>& MaterialParam::epsZ() { return _eps_z; }

xt::xarray<double>& MaterialParam::muX() { return _mu_x; }

xt::xarray<double>& MaterialParam::muY() { return _mu_y; }

xt::xarray<double>& MaterialParam::muZ() { return _mu_z; }

xt::xarray<double>& MaterialParam::sigmaEX() { return _sigma_e_x; }

xt::xarray<double>& MaterialParam::sigmaEY() { return _sigma_e_y; }

xt::xarray<double>& MaterialParam::sigmaEZ() { return _sigma_e_z; }

xt::xarray<double>& MaterialParam::sigmaMX() { return _sigma_m_x; }

xt::xarray<double>& MaterialParam::sigmaMY() { return _sigma_m_y; }

xt::xarray<double>& MaterialParam::sigmaMZ() { return _sigma_m_z; }

xt::xarray<double>& MaterialParam::property(MaterialParam::Property property) {
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

xt::xarray<double>& MaterialParam::property(MaterialParam::Attribute attribute,
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

void MaterialParam::allocate(std::size_t nx, std::size_t ny, std::size_t nz) {
  _eps_x = xt::zeros<double>({nx, ny + 1, nz + 1});
  _eps_y = xt::zeros<double>({nx + 1, ny, nz + 1});
  _eps_z = xt::zeros<double>({nx + 1, ny + 1, nz});
  _eps_x.fill(constant::EPSILON_0);
  _eps_y.fill(constant::EPSILON_0);
  _eps_z.fill(constant::EPSILON_0);

  _mu_x = xt::zeros<double>({nx + 1, ny, nz});
  _mu_y = xt::zeros<double>({nx, ny + 1, nz});
  _mu_z = xt::zeros<double>({nx, ny, nz + 1});
  _mu_x.fill(constant::MU_0);
  _mu_y.fill(constant::MU_0);
  _mu_z.fill(constant::MU_0);

  _sigma_e_x = xt::zeros<double>({nx, ny + 1, nz + 1});
  _sigma_e_y = xt::zeros<double>({nx + 1, ny, nz + 1});
  _sigma_e_z = xt::zeros<double>({nx + 1, ny + 1, nz});
  _sigma_e_x.fill(constant::SIGMA_E_ZERO_APPROX);
  _sigma_e_y.fill(constant::SIGMA_E_ZERO_APPROX);
  _sigma_e_z.fill(constant::SIGMA_E_ZERO_APPROX);

  _sigma_m_x = xt::zeros<double>({nx + 1, ny, nz});
  _sigma_m_y = xt::zeros<double>({nx, ny + 1, nz});
  _sigma_m_z = xt::zeros<double>({nx, ny, nz + 1});
  _sigma_m_x.fill(constant::SIGMA_M_ZERO_APPROX);
  _sigma_m_y.fill(constant::SIGMA_M_ZERO_APPROX);
  _sigma_m_z.fill(constant::SIGMA_M_ZERO_APPROX);
}

}  // namespace xfdtd
