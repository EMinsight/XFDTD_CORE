#ifndef _XFDTD_CORE_CALCULATION_PARAM_H_
#define _XFDTD_CORE_CALCULATION_PARAM_H_

#include <xfdtd/calculation_param/time_param.h>
#include <xfdtd/common/type_define.h>
#include <xfdtd/coordinate_system/coordinate_system.h>
#include <xfdtd/electromagnetic_field/electromagnetic_field.h>
#include <xfdtd/exception/exception.h>
#include <xfdtd/util/transform.h>

#include <memory>
#include <vector>

namespace xfdtd {

class GridSpace;

// forward declare
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

class FDTDUpdateCoefficient {
 public:
  template <EMF::Attribute c, Axis::XYZ xyz>
  auto coeff() const -> const Array3D<Real>&;

  template <EMF::Attribute c, Axis::XYZ xyz>
  auto coeff() -> Array3D<Real>&;

  template <EMF::Attribute a, Axis::XYZ xyz_a, EMF::Attribute b,
            Axis::XYZ xyz_b>
  auto coeff() const -> const Array3D<Real>&;

  template <EMF::Attribute a, Axis::XYZ xyz_a, EMF::Attribute b,
            Axis::XYZ xyz_b>
  auto coeff() -> Array3D<Real>&;

  template <EMF::Attribute c, Axis::XYZ xyz>
  auto coeffDualA() const -> const Array3D<Real>&;

  template <EMF::Attribute c, Axis::XYZ xyz>
  auto coeffDualA() -> Array3D<Real>&;

  template <EMF::Attribute c, Axis::XYZ xyz>
  auto coeffDualB() const -> const Array3D<Real>&;

  template <EMF::Attribute c, Axis::XYZ xyz>
  auto coeffDualB() -> Array3D<Real>&;

  auto save(const std::string& dir) const -> void;

 public:
  const Array3D<Real>& cexe() const;

  const Array3D<Real>& cexhy() const;

  const Array3D<Real>& cexhz() const;

  const Array3D<Real>& ceye() const;

  const Array3D<Real>& ceyhz() const;

  const Array3D<Real>& ceyhx() const;

  const Array3D<Real>& ceze() const;

  const Array3D<Real>& cezhx() const;

  const Array3D<Real>& cezhy() const;

  const Array3D<Real>& chxh() const;

  const Array3D<Real>& chxey() const;

  const Array3D<Real>& chxez() const;

  const Array3D<Real>& chyh() const;

  const Array3D<Real>& chyez() const;

  const Array3D<Real>& chyex() const;

  const Array3D<Real>& chzh() const;

  const Array3D<Real>& chzex() const;

  const Array3D<Real>& chzey() const;

  Array3D<Real>& cexe();

  Array3D<Real>& cexhy();

  Array3D<Real>& cexhz();

  Array3D<Real>& ceye();

  Array3D<Real>& ceyhz();

  Array3D<Real>& ceyhx();

  Array3D<Real>& ceze();

  Array3D<Real>& cezhx();

  Array3D<Real>& cezhy();

  Array3D<Real>& chxh();

  Array3D<Real>& chxey();

  Array3D<Real>& chxez();

  Array3D<Real>& chyh();

  Array3D<Real>& chyez();

  Array3D<Real>& chyex();

  Array3D<Real>& chzh();

  Array3D<Real>& chzex();

  Array3D<Real>& chzey();

 private:
  Array3D<Real> _cexe;
  Array3D<Real> _cexhy;
  Array3D<Real> _cexhz;
  Array3D<Real> _ceye;
  Array3D<Real> _ceyhz;
  Array3D<Real> _ceyhx;
  Array3D<Real> _ceze;
  Array3D<Real> _cezhx;
  Array3D<Real> _cezhy;

  Array3D<Real> _chxh;
  Array3D<Real> _chxey;
  Array3D<Real> _chxez;
  Array3D<Real> _chyh;
  Array3D<Real> _chyez;
  Array3D<Real> _chyex;
  Array3D<Real> _chzh;
  Array3D<Real> _chzex;
  Array3D<Real> _chzey;
};

template <EMF::Attribute c, Axis::XYZ xyz>
inline auto FDTDUpdateCoefficient::coeff() const -> const Array3D<Real>& {
  constexpr auto f = transform::attributeXYZToField<c, xyz>();
  switch (f) {
    case EMF::Field::EX:
      return cexe();
    case EMF::Field::EY:
      return ceye();
    case EMF::Field::EZ:
      return ceze();
    case EMF::Field::HX:
      return chxh();
    case EMF::Field::HY:
      return chyh();
    case EMF::Field::HZ:
      return chzh();
    default:
      throw XFDTDException{
          "FDTDUpdateCoefficient::coeff(): Invalid EMF::Field"};
  }
}

template <EMF::Attribute c, Axis::XYZ xyz>
inline auto FDTDUpdateCoefficient::coeff() -> Array3D<Real>& {
  return const_cast<Array3D<Real>&>(
      static_cast<const FDTDUpdateCoefficient*>(this)->coeff<c, xyz>());
}

template <EMF::Attribute a, Axis::XYZ xyz_a, EMF::Attribute b, Axis::XYZ xyz_b>
auto FDTDUpdateCoefficient::coeff() const -> const Array3D<Real>& {
  constexpr auto f = transform::attributeXYZToField<a, xyz_a>();
  constexpr auto g = transform::attributeXYZToField<b, xyz_b>();
  if constexpr (f == EMF::Field::EX && g == EMF::Field::HY) {
    return cexhy();
  } else if constexpr (f == EMF::Field::EX && g == EMF::Field::HZ) {
    return cexhz();
  } else if constexpr (f == EMF::Field::EY && g == EMF::Field::HZ) {
    return ceyhz();
  } else if constexpr (f == EMF::Field::EY && g == EMF::Field::HX) {
    return ceyhx();
  } else if constexpr (f == EMF::Field::EZ && g == EMF::Field::HX) {
    return cezhx();
  } else if constexpr (f == EMF::Field::EZ && g == EMF::Field::HY) {
    return cezhy();
  } else if constexpr (f == EMF::Field::HX && g == EMF::Field::EY) {
    return chxey();
  } else if constexpr (f == EMF::Field::HX && g == EMF::Field::EZ) {
    return chxez();
  } else if constexpr (f == EMF::Field::HY && g == EMF::Field::EZ) {
    return chyez();
  } else if constexpr (f == EMF::Field::HY && g == EMF::Field::EX) {
    return chyex();
  } else if constexpr (f == EMF::Field::HZ && g == EMF::Field::EX) {
    return chzex();
  } else if constexpr (f == EMF::Field::HZ && g == EMF::Field::EY) {
    return chzey();
  } else {
    throw XFDTDException{"FDTDUpdateCoefficient::coeff(): Invalid EMF::Field"};
  }
}

template <EMF::Attribute a, Axis::XYZ xyz_a, EMF::Attribute b, Axis::XYZ xyz_b>
auto FDTDUpdateCoefficient::coeff() -> Array3D<Real>& {
  return const_cast<Array3D<Real>&>(
      static_cast<const FDTDUpdateCoefficient*>(this)
          ->coeff<a, xyz_a, b, xyz_b>());
}

template <EMF::Attribute c, Axis::XYZ xyz>
inline auto FDTDUpdateCoefficient::coeffDualA() const -> const Array3D<Real>& {
  constexpr auto a = EMF::dualAttribute(c);
  const auto xyz_a = Axis::tangentialAAxis<xyz>();
  return coeff<c, xyz, a, xyz_a>();
}

template <EMF::Attribute c, Axis::XYZ xyz>
inline auto FDTDUpdateCoefficient::coeffDualA() -> Array3D<Real>& {
  return const_cast<Array3D<Real>&>(
      static_cast<const FDTDUpdateCoefficient*>(this)->coeffDualA<c, xyz>());
}

template <EMF::Attribute c, Axis::XYZ xyz>
inline auto FDTDUpdateCoefficient::coeffDualB() const -> const Array3D<Real>& {
  constexpr auto a = EMF::dualAttribute(c);
  const auto xyz_b = Axis::tangentialBAxis<xyz>();
  return coeff<c, xyz, a, xyz_b>();
}

template <EMF::Attribute c, Axis::XYZ xyz>
inline auto FDTDUpdateCoefficient::coeffDualB() -> Array3D<Real>& {
  return const_cast<Array3D<Real>&>(
      static_cast<const FDTDUpdateCoefficient*>(this)->coeffDualB<c, xyz>());
}

class XFDTDCalculationParamException : public XFDTDException {
 public:
  explicit XFDTDCalculationParamException(std::string message);

  ~XFDTDCalculationParamException() override = default;
};

class CalculationParam {
 public:
  CalculationParam();

  ~CalculationParam() = default;

  const std::unique_ptr<TimeParam>& timeParam() const;

  const std::unique_ptr<MaterialParam>& materialParam() const;

  const std::unique_ptr<FDTDUpdateCoefficient>& fdtdCoefficient() const;

  std::unique_ptr<TimeParam>& timeParam();

  std::unique_ptr<MaterialParam>& materialParam();

  std::unique_ptr<FDTDUpdateCoefficient>& fdtdCoefficient();

  void generateMaterialSpaceParam(const GridSpace* grid_space);

  void calculateCoefficient(const GridSpace* grid_space);

  void setMaterialParam(std::unique_ptr<MaterialParam> material_param);

  void setTimeParam(std::unique_ptr<TimeParam> time_param);

 private:
  std::unique_ptr<TimeParam> _time_param;
  std::unique_ptr<MaterialParam> _material_param;
  std::unique_ptr<FDTDUpdateCoefficient> _fdtd_coefficient;
};

}  // namespace xfdtd

#endif  // _XFDTD_CORE_CALCULATION_PARAM_H_
