#ifndef __XFDTD_CORE_FDTD_UPDATE_COEFFICIENT_H__
#define __XFDTD_CORE_FDTD_UPDATE_COEFFICIENT_H__

#include "xfdtd/common/type_define.h"
#include "xfdtd/coordinate_system/coordinate_system.h"
#include "xfdtd/electromagnetic_field/electromagnetic_field.h"
#include "xfdtd/util/transform.h"

namespace xfdtd {

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
  if constexpr (f == EMF::Field::EX) {
    return cexe();
  } else if constexpr (f == EMF::Field::EY) {
    return ceye();
  } else if constexpr (f == EMF::Field::EZ) {
    return ceze();
  } else if constexpr (f == EMF::Field::HX) {
    return chxh();
  } else if constexpr (f == EMF::Field::HY) {
    return chyh();
  } else if constexpr (f == EMF::Field::HZ) {
    return chzh();
  } else {
    static_assert(f == EMF::Field::EX || f == EMF::Field::EY ||
                      f == EMF::Field::EZ || f == EMF::Field::HX ||
                      f == EMF::Field::HY || f == EMF::Field::HZ,
                  "FDTDUpdateCoefficient::coeff(): Invalid EMF::Field");
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
    static_assert(f == EMF::Field::EX || f == EMF::Field::EY ||
                      f == EMF::Field::EZ || f == EMF::Field::HX ||
                      f == EMF::Field::HY || f == EMF::Field::HZ,
                  "FDTDUpdateCoefficient::coeff(): Invalid EMF::Field");
  }
}

template <EMF::Attribute a, Axis::XYZ xyz_a, EMF::Attribute b, Axis::XYZ xyz_b>
auto FDTDUpdateCoefficient::coeff() -> Array3D<Real>& {
  return const_cast<Array3D<Real>&>(
      static_cast<const FDTDUpdateCoefficient*>(this)
          ->coeff<a, xyz_a, b, xyz_b>());
}

}  // namespace xfdtd

#endif  // __XFDTD_CORE_FDTD_UPDATE_COEFFICIENT_H__
