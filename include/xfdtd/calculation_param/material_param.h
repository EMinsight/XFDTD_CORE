#ifndef __XFDTD_CORE_MATERIAL_PARAM_H__
#define __XFDTD_CORE_MATERIAL_PARAM_H__

#include <xfdtd/common/type_define.h>
#include <xfdtd/coordinate_system/coordinate_system.h>

#include <vector>

namespace xfdtd {

class Material;

class MaterialParam {
 public:
  enum class Attribute { EPSILON, MU, SIGMA_E, SIGMA_M };

  enum class Property {
    EPS_X,
    EPS_Y,
    EPS_Z,
    MU_X,
    MU_Y,
    MU_Z,
    SIGMA_E_X,
    SIGMA_E_Y,
    SIGMA_E_Z,
    SIGMA_M_X,
    SIGMA_M_Y,
    SIGMA_M_Z
  };

 public:
  static Attribute fromPropertyToAttribute(Property property);

  static Axis::XYZ fromPropertyToXYZ(Property property);

  static Property fromAttributeAndXYZ(Attribute attribute, Axis::XYZ xyz);

  const Array3D<Real>& epsX() const;

  const Array3D<Real>& epsY() const;

  const Array3D<Real>& epsZ() const;

  const Array3D<Real>& muX() const;

  const Array3D<Real>& muY() const;

  const Array3D<Real>& muZ() const;

  const Array3D<Real>& sigmaEX() const;

  const Array3D<Real>& sigmaEY() const;

  const Array3D<Real>& sigmaEZ() const;

  const Array3D<Real>& sigmaMX() const;

  const Array3D<Real>& sigmaMY() const;

  const Array3D<Real>& sigmaMZ() const;

  const Array3D<Real>& property(Property property) const;

  const Array3D<Real>& property(Attribute attribute, Axis::XYZ xyz) const;

  const auto& materialArray() const;

  Array3D<Real>& epsX();

  Array3D<Real>& epsY();

  Array3D<Real>& epsZ();

  Array3D<Real>& muX();

  Array3D<Real>& muY();

  Array3D<Real>& muZ();

  Array3D<Real>& sigmaEX();

  Array3D<Real>& sigmaEY();

  Array3D<Real>& sigmaEZ();

  Array3D<Real>& sigmaMX();

  Array3D<Real>& sigmaMY();

  Array3D<Real>& sigmaMZ();

  Array3D<Real>& property(Property property);

  Array3D<Real>& property(Attribute attribute, Axis::XYZ xyz);

  auto&& materialArray();

  void addMaterial(std::shared_ptr<Material> material);

  void allocate(std::size_t nx, std::size_t ny, std::size_t nz);

 private:
  Array3D<Real> _eps_x;
  Array3D<Real> _eps_y;
  Array3D<Real> _eps_z;
  Array3D<Real> _mu_x;
  Array3D<Real> _mu_y;
  Array3D<Real> _mu_z;
  Array3D<Real> _sigma_e_x;
  Array3D<Real> _sigma_e_y;
  Array3D<Real> _sigma_e_z;
  Array3D<Real> _sigma_m_x;
  Array3D<Real> _sigma_m_y;
  Array3D<Real> _sigma_m_z;

  std::vector<std::shared_ptr<Material>> _materials;
};

inline const auto& MaterialParam::materialArray() const { return _materials; }

inline auto&& MaterialParam::materialArray() { return _materials; }

}  // namespace xfdtd

#endif  // __XFDTD_CORE_MATERIAL_PARAM_H__
