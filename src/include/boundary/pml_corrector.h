#ifndef _XFDTD_CORE_PML_CORRECTOR_H_
#define _XFDTD_CORE_PML_CORRECTOR_H_

#include <xfdtd/common/type_define.h>

#include <cstddef>
#include <xtensor/xarray.hpp>

#include "corrector/corrector.h"

namespace xfdtd {

class PMLCorrector : public Corrector {
 public:
  PMLCorrector(std::size_t pml_global_e_start, std::size_t pml_global_h_start,
               std::size_t pml_node_e_start, std::size_t pml_node_h_start,
               std::size_t a_s, std::size_t a_n, std::size_t b_s,
               std::size_t b_n, std::size_t c_e_s, std::size_t c_e_n,
               std::size_t c_h_s, std::size_t c_h_n, std::size_t offset_c,
               Array1D<Real>& coeff_a_e, Array1D<Real>& coeff_b_e,
               Array1D<Real>& coeff_a_h, Array1D<Real>& coeff_b_h,
               Array3D<Real>& c_ea_psi_hb, Array3D<Real>& c_eb_psi_ha,
               Array3D<Real>& c_ha_psi_eb, Array3D<Real>& c_hb_psi_ea,
               Array3D<Real>& ea_psi_hb, Array3D<Real>& eb_psi_ha,
               Array3D<Real>& ha_psi_eb, Array3D<Real>& hb_psi_ea,
               Array3D<Real>& ea, Array3D<Real>& eb, Array3D<Real>& ha,
               Array3D<Real>& hb)
      : _pml_global_e_start{pml_global_e_start},
        _pml_global_h_start{pml_global_h_start},
        _pml_node_e_start{pml_node_e_start},
        _pml_node_h_start{pml_node_h_start},
        _a_s{a_s},
        _a_n{a_n},
        _b_s{b_s},
        _b_n{b_n},
        _c_e_s{c_e_s},
        _c_e_n{c_e_n},
        _c_h_s{c_h_s},
        _c_h_n{c_h_n},
        _offset_c{offset_c},
        _coeff_a_e{coeff_a_e},
        _coeff_b_e{coeff_b_e},
        _coeff_a_h{coeff_a_h},
        _coeff_b_h{coeff_b_h},
        _c_ea_psi_hb{c_ea_psi_hb},
        _c_eb_psi_ha{c_eb_psi_ha},
        _c_ha_psi_eb{c_ha_psi_eb},
        _c_hb_psi_ea{c_hb_psi_ea},
        _ea_psi_hb{ea_psi_hb},
        _eb_psi_ha{eb_psi_ha},
        _ha_psi_eb{ha_psi_eb},
        _hb_psi_ea{hb_psi_ea},
        _ea{ea},
        _eb{eb},
        _ha{ha},
        _hb{hb} {}

  ~PMLCorrector() override = default;

 protected:
  std::size_t _pml_global_e_start, _pml_global_h_start;
  std::size_t _pml_node_e_start, _pml_node_h_start;
  std::size_t _offset_c;
  std::size_t _a_s, _a_n;
  std::size_t _b_s, _b_n;
  std::size_t _c_e_s, _c_e_n;
  std::size_t _c_h_s, _c_h_n;
  Array1D<Real>&_coeff_a_e, &_coeff_b_e, &_coeff_a_h, &_coeff_b_h;
  Array3D<Real>&_c_ea_psi_hb, &_c_eb_psi_ha, &_c_ha_psi_eb, &_c_hb_psi_ea;
  Array3D<Real>&_ea_psi_hb, &_eb_psi_ha, &_ha_psi_eb, &_hb_psi_ea;
  Array3D<Real>&_ea, &_eb, &_ha, &_hb;
};

class PMLCorrectorX : public PMLCorrector {
 public:
  PMLCorrectorX(std::size_t pml_global_e_start, std::size_t pml_global_h_start,
                std::size_t pml_node_e_start, std::size_t pml_node_h_start,
                std::size_t a_s, std::size_t a_n, std::size_t b_s,
                std::size_t b_n, std::size_t c_e_s, std::size_t c_e_n,
                std::size_t c_h_s, std::size_t c_h_n, std::size_t offset_c,
                Array1D<Real>& coeff_a_e, Array1D<Real>& coeff_b_e,
                Array1D<Real>& coeff_a_h, Array1D<Real>& coeff_b_h,
                Array3D<Real>& c_ea_psi_hb, Array3D<Real>& c_eb_psi_ha,
                Array3D<Real>& c_ha_psi_eb, Array3D<Real>& c_hb_psi_ea,
                Array3D<Real>& ea_psi_hb, Array3D<Real>& eb_psi_ha,
                Array3D<Real>& ha_psi_eb, Array3D<Real>& hb_psi_ea,
                Array3D<Real>& ea, Array3D<Real>& eb, Array3D<Real>& ha,
                Array3D<Real>& hb)
      : PMLCorrector(pml_global_e_start, pml_global_h_start, pml_node_e_start,
                     pml_node_h_start, a_s, a_n, b_s, b_n, c_e_s, c_e_n, c_h_s,
                     c_h_n, offset_c, coeff_a_e, coeff_b_e, coeff_a_h,
                     coeff_b_h, c_ea_psi_hb, c_eb_psi_ha, c_ha_psi_eb,
                     c_hb_psi_ea, ea_psi_hb, eb_psi_ha, ha_psi_eb, hb_psi_ea,
                     ea, eb, ha, hb) {}

  ~PMLCorrectorX() override = default;

  void correctE() override;

  void correctH() override;

  std::string toString() const override;

 private:
};

class PMLCorrectorY : public PMLCorrector {
 public:
  PMLCorrectorY(std::size_t pml_global_e_start, std::size_t pml_global_h_start,
                std::size_t pml_node_e_start, std::size_t pml_node_h_start,
                std::size_t a_s, std::size_t a_n, std::size_t b_s,
                std::size_t b_n, std::size_t c_e_s, std::size_t c_e_n,
                std::size_t c_h_s, std::size_t c_h_n, std::size_t offset_c,
                Array1D<Real>& coeff_a_e, Array1D<Real>& coeff_b_e,
                Array1D<Real>& coeff_a_h, Array1D<Real>& coeff_b_h,
                Array3D<Real>& c_ea_psi_hb, Array3D<Real>& c_eb_psi_ha,
                Array3D<Real>& c_ha_psi_eb, Array3D<Real>& c_hb_psi_ea,
                Array3D<Real>& ea_psi_hb, Array3D<Real>& eb_psi_ha,
                Array3D<Real>& ha_psi_eb, Array3D<Real>& hb_psi_ea,
                Array3D<Real>& ea, Array3D<Real>& eb, Array3D<Real>& ha,
                Array3D<Real>& hb)
      : PMLCorrector(pml_global_e_start, pml_global_h_start, pml_node_e_start,
                     pml_node_h_start, a_s, a_n, b_s, b_n, c_e_s, c_e_n, c_h_s,
                     c_h_n, offset_c, coeff_a_e, coeff_b_e, coeff_a_h,
                     coeff_b_h, c_ea_psi_hb, c_eb_psi_ha, c_ha_psi_eb,
                     c_hb_psi_ea, ea_psi_hb, eb_psi_ha, ha_psi_eb, hb_psi_ea,
                     ea, eb, ha, hb) {}

  ~PMLCorrectorY() override = default;

  void correctE() override;

  void correctH() override;

  std::string toString() const override;

 private:
};

class PMLCorrectorZ : public PMLCorrector {
 public:
  PMLCorrectorZ(std::size_t pml_global_e_start, std::size_t pml_global_h_start,
                std::size_t pml_node_e_start, std::size_t pml_node_h_start,
                std::size_t a_s, std::size_t a_n, std::size_t b_s,
                std::size_t b_n, std::size_t c_e_s, std::size_t c_e_n,
                std::size_t c_h_s, std::size_t c_h_n, std::size_t offset_c,
                Array1D<Real>& coeff_a_e, Array1D<Real>& coeff_b_e,
                Array1D<Real>& coeff_a_h, Array1D<Real>& coeff_b_h,
                Array3D<Real>& c_ea_psi_hb, Array3D<Real>& c_eb_psi_ha,
                Array3D<Real>& c_ha_psi_eb, Array3D<Real>& c_hb_psi_ea,
                Array3D<Real>& ea_psi_hb, Array3D<Real>& eb_psi_ha,
                Array3D<Real>& ha_psi_eb, Array3D<Real>& hb_psi_ea,
                Array3D<Real>& ea, Array3D<Real>& eb, Array3D<Real>& ha,
                Array3D<Real>& hb)
      : PMLCorrector(pml_global_e_start, pml_global_h_start, pml_node_e_start,
                     pml_node_h_start, a_s, a_n, b_s, b_n, c_e_s, c_e_n, c_h_s,
                     c_h_n, offset_c, coeff_a_e, coeff_b_e, coeff_a_h,
                     coeff_b_h, c_ea_psi_hb, c_eb_psi_ha, c_ha_psi_eb,
                     c_hb_psi_ea, ea_psi_hb, eb_psi_ha, ha_psi_eb, hb_psi_ea,
                     ea, eb, ha, hb) {}

  ~PMLCorrectorZ() override = default;

  void correctE() override;

  void correctH() override;

  std::string toString() const override;

 private:
};

}  // namespace xfdtd

#endif  // _XFDTD_CORE_PML_CORRECTOR_H_
