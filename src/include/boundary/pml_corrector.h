#ifndef _XFDTD_LIB_PML_CORRECTOR_H_
#define _XFDTD_LIB_PML_CORRECTOR_H_

#include <cstddef>
#include <xtensor/xarray.hpp>

#include "corrector/corrector.h"

namespace xfdtd {

class PMLCorrector : public Corrector {
 public:
  PMLCorrector(std::size_t pml_e_c_start, std::size_t pml_h_c_start,
               std::size_t a_s, std::size_t a_n, std::size_t b_s,
               std::size_t b_n, std::size_t c_e_s, std::size_t c_e_n,
               std::size_t c_h_s, std::size_t c_h_n,
               xt::xarray<double>& coeff_a_e, xt::xarray<double>& coeff_b_e,
               xt::xarray<double>& coeff_a_h, xt::xarray<double>& coeff_b_h,
               xt::xarray<double>& c_ea_psi_hb, xt::xarray<double>& c_eb_psi_ha,
               xt::xarray<double>& c_ha_psi_eb, xt::xarray<double>& c_hb_psi_ea,
               xt::xarray<double>& ea_psi_hb, xt::xarray<double>& eb_psi_ha,
               xt::xarray<double>& ha_psi_eb, xt::xarray<double>& hb_psi_ea,
               xt::xarray<double>& ea, xt::xarray<double>& eb,
               xt::xarray<double>& ha, xt::xarray<double>& hb)
      : _pml_e_c_start{pml_e_c_start},
        _pml_h_c_start{pml_h_c_start},
        _a_s{a_s},
        _a_n{a_n},
        _b_s{b_s},
        _b_n{b_n},
        _c_e_s{c_e_s},
        _c_e_n{c_e_n},
        _c_h_s{c_h_s},
        _c_h_n{c_h_n},
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
  std::size_t _pml_e_c_start, _pml_h_c_start;
  std::size_t _a_s, _a_n;
  std::size_t _b_s, _b_n;
  std::size_t _c_e_s, _c_e_n;
  std::size_t _c_h_s, _c_h_n;
  xt::xarray<double>&_coeff_a_e, &_coeff_b_e, &_coeff_a_h, &_coeff_b_h;
  xt::xarray<double>&_c_ea_psi_hb, &_c_eb_psi_ha, &_c_ha_psi_eb, &_c_hb_psi_ea;
  xt::xarray<double>&_ea_psi_hb, &_eb_psi_ha, &_ha_psi_eb, &_hb_psi_ea;
  xt::xarray<double>&_ea, &_eb, &_ha, &_hb;
};

class PMLCorrectorX : public PMLCorrector {
 public:
  PMLCorrectorX(std::size_t pml_e_c_start, std::size_t pml_h_c_start,
                std::size_t a_s, std::size_t a_n, std::size_t b_s,
                std::size_t b_n, std::size_t c_e_s, std::size_t c_e_n,
                std::size_t c_h_s, std::size_t c_h_n,
                xt::xarray<double>& coeff_a_e, xt::xarray<double>& coeff_b_e,
                xt::xarray<double>& coeff_a_h, xt::xarray<double>& coeff_b_h,
                xt::xarray<double>& c_ea_psi_hb,
                xt::xarray<double>& c_eb_psi_ha,
                xt::xarray<double>& c_ha_psi_eb,
                xt::xarray<double>& c_hb_psi_ea, xt::xarray<double>& ea_psi_hb,
                xt::xarray<double>& eb_psi_ha, xt::xarray<double>& ha_psi_eb,
                xt::xarray<double>& hb_psi_ea, xt::xarray<double>& ea,
                xt::xarray<double>& eb, xt::xarray<double>& ha,
                xt::xarray<double>& hb)
      : PMLCorrector(pml_e_c_start, pml_h_c_start, a_s, a_n, b_s, b_n, c_e_s,
                     c_e_n, c_h_s, c_h_n, coeff_a_e, coeff_b_e, coeff_a_h,
                     coeff_b_h, c_ea_psi_hb, c_eb_psi_ha, c_ha_psi_eb,
                     c_hb_psi_ea, ea_psi_hb, eb_psi_ha, ha_psi_eb, hb_psi_ea,
                     ea, eb, ha, hb) {}

  ~PMLCorrectorX() override = default;

  void correctE() override;

  void correctH() override;

 private:
};

class PMLCorrectorY : public PMLCorrector {
 public:
  PMLCorrectorY(std::size_t pml_e_c_start, std::size_t pml_h_c_start,
                std::size_t a_s, std::size_t a_n, std::size_t b_s,
                std::size_t b_n, std::size_t c_e_s, std::size_t c_e_n,
                std::size_t c_h_s, std::size_t c_h_n,
                xt::xarray<double>& coeff_a_e, xt::xarray<double>& coeff_b_e,
                xt::xarray<double>& coeff_a_h, xt::xarray<double>& coeff_b_h,
                xt::xarray<double>& c_ea_psi_hb,
                xt::xarray<double>& c_eb_psi_ha,
                xt::xarray<double>& c_ha_psi_eb,
                xt::xarray<double>& c_hb_psi_ea, xt::xarray<double>& ea_psi_hb,
                xt::xarray<double>& eb_psi_ha, xt::xarray<double>& ha_psi_eb,
                xt::xarray<double>& hb_psi_ea, xt::xarray<double>& ea,
                xt::xarray<double>& eb, xt::xarray<double>& ha,
                xt::xarray<double>& hb)
      : PMLCorrector(pml_e_c_start, pml_h_c_start, a_s, a_n, b_s, b_n, c_e_s,
                     c_e_n, c_h_s, c_h_n, coeff_a_e, coeff_b_e, coeff_a_h,
                     coeff_b_h, c_ea_psi_hb, c_eb_psi_ha, c_ha_psi_eb,
                     c_hb_psi_ea, ea_psi_hb, eb_psi_ha, ha_psi_eb, hb_psi_ea,
                     ea, eb, ha, hb) {}

  ~PMLCorrectorY() override = default;

  void correctE() override;

  void correctH() override;

 private:
};

class PMLCorrectorZ : public PMLCorrector {
 public:
  PMLCorrectorZ(std::size_t pml_e_c_start, std::size_t pml_h_c_start,
                std::size_t a_s, std::size_t a_n, std::size_t b_s,
                std::size_t b_n, std::size_t c_e_s, std::size_t c_e_n,
                std::size_t c_h_s, std::size_t c_h_n,
                xt::xarray<double>& coeff_a_e, xt::xarray<double>& coeff_b_e,
                xt::xarray<double>& coeff_a_h, xt::xarray<double>& coeff_b_h,
                xt::xarray<double>& c_ea_psi_hb,
                xt::xarray<double>& c_eb_psi_ha,
                xt::xarray<double>& c_ha_psi_eb,
                xt::xarray<double>& c_hb_psi_ea, xt::xarray<double>& ea_psi_hb,
                xt::xarray<double>& eb_psi_ha, xt::xarray<double>& ha_psi_eb,
                xt::xarray<double>& hb_psi_ea, xt::xarray<double>& ea,
                xt::xarray<double>& eb, xt::xarray<double>& ha,
                xt::xarray<double>& hb)
      : PMLCorrector(pml_e_c_start, pml_h_c_start, a_s, a_n, b_s, b_n, c_e_s,
                     c_e_n, c_h_s, c_h_n, coeff_a_e, coeff_b_e, coeff_a_h,
                     coeff_b_h, c_ea_psi_hb, c_eb_psi_ha, c_ha_psi_eb,
                     c_hb_psi_ea, ea_psi_hb, eb_psi_ha, ha_psi_eb, hb_psi_ea,
                     ea, eb, ha, hb) {}

  ~PMLCorrectorZ() override = default;

  void correctE() override;

  void correctH() override;

 private:
};

}  // namespace xfdtd

#endif  // _XFDTD_LIB_PML_CORRECTOR_H_
