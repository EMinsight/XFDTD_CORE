#ifndef __XFDTD_CORE_FDTD_UPDATE_COEFFICIENT_H__
#define __XFDTD_CORE_FDTD_UPDATE_COEFFICIENT_H__

#include "xfdtd/common/type_define.h"
#include "xfdtd/coordinate_system/coordinate_system.h"
#include "xfdtd/electromagnetic_field/electromagnetic_field.h"

namespace xfdtd {

class FDTDUpdateCoefficient {
 public:
  template <EMF::Attribute attribute, Axis::XYZ>
  auto coefficent() const {
    if constexpr (attribute == EMF::Attribute::E) {
      if constexpr (xyz == Axis::XYZ::X) {
        return cexe();
      } else if constexpr (xyz == Axis::XYZ::Y) {
        return ceye();
      } else if constexpr (xyz == Axis::XYZ::Z) {
        return ceze();
      }
    } else {
      if constexpr (xyz == Axis::XYZ::X) {
        return chxh();
      } else if constexpr (xyz == Axis::XYZ::Y) {
        return chyh();
      } else if constexpr (xyz == Axis::XYZ::Z) {
        return chzh();
      }
    }
  }

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

}  // namespace xfdtd

#endif  // __XFDTD_CORE_FDTD_UPDATE_COEFFICIENT_H__
