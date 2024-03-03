#include "boundary/pml_corrector.h"

#include "boundary/pml_scheme.h"

namespace xfdtd {

void PMLCorrectorX::correctE() {
  for (std::size_t i{_c_e_s}; i < _c_e_s + _c_e_n; ++i) {
    auto ii{i - _pml_e_c_start};
    for (std::size_t j{_a_s}; j < _a_s + _a_n; ++j) {
      for (std::size_t k{_b_s}; k < _b_s + _b_n + 1; ++k) {
        correctPML(_ea(i, j, k), _ea_psi_hb(ii, j, k), _coeff_a_e(ii),
                   _coeff_b_e(ii), _hb(i, j, k), _hb(i - 1, j, k),
                   _c_ea_psi_hb(ii, j, k));
      }
    }

    for (std::size_t j{_a_s}; j < _a_s + _a_n + 1; ++j) {
      for (std::size_t k{_b_s}; k < _b_s + _b_n; ++k) {
        correctPML(_eb(i, j, k), _eb_psi_ha(ii, j, k), _coeff_a_e(ii),
                   _coeff_b_e(ii), _ha(i, j, k), _ha(i - 1, j, k),
                   _c_eb_psi_ha(ii, j, k));
      }
    }
  }
}

void PMLCorrectorX::correctH() {
  for (std::size_t i{_c_h_s}; i < _c_h_s + _c_h_n; ++i) {
    auto ii{i - _pml_h_c_start};
    for (std::size_t j{_a_s}; j < _a_s + _a_n + 1; ++j) {
      for (std::size_t k{_b_s}; k < _b_s + _b_n; ++k) {
        correctPML(_ha(i, j, k), _ha_psi_eb(ii, j, k), _coeff_a_h(ii),
                   _coeff_b_h(ii), _eb(i + 1, j, k), _eb(i, j, k),
                   _c_ha_psi_eb(ii, j, k));
      }
    }

    for (std::size_t j{_a_s}; j < _a_s + _a_n; ++j) {
      for (std::size_t k{_b_s}; k < _b_s + _b_n + 1; ++k) {
        correctPML(_hb(i, j, k), _hb_psi_ea(ii, j, k), _coeff_a_h(ii),
                   _coeff_b_h(ii), _ea(i + 1, j, k), _ea(i, j, k),
                   _c_hb_psi_ea(ii, j, k));
      }
    }
  }
}

void PMLCorrectorY::correctE() {
  for (std::size_t i{_b_s}; i < _b_s + _b_n + 1; ++i) {
    for (std::size_t j{_c_e_s}; j < _c_e_s + _c_e_n; ++j) {
      auto jj{j - _pml_e_c_start};
      for (std::size_t k{_a_s}; k < _a_s + _a_n; ++k) {
        correctPML(_ea(i, j, k), _ea_psi_hb(i, jj, k), _coeff_a_e(jj),
                   _coeff_b_e(jj), _hb(i, j, k), _hb(i, j - 1, k),
                   _c_ea_psi_hb(i, jj, k));
      }
    }
  }

  for (std::size_t i{_b_s}; i < _b_s + _b_n; ++i) {
    for (std::size_t j{_c_e_s}; j < _c_e_s + _c_e_n; ++j) {
      auto jj{j - _pml_e_c_start};
      for (std::size_t k{_a_s}; k < _a_s + _a_n + 1; ++k) {
        correctPML(_eb(i, j, k), _eb_psi_ha(i, jj, k), _coeff_a_e(jj),
                   _coeff_b_e(jj), _ha(i, j, k), _ha(i, j - 1, k),
                   _c_eb_psi_ha(i, jj, k));
      }
    }
  }
}

void PMLCorrectorY::correctH() {
  for (std::size_t i{_b_s}; i < _b_s + _b_n; ++i) {
    for (std::size_t j{_c_h_s}; j < _c_h_s + _c_h_n; ++j) {
      for (std::size_t k{_a_s}; k < _a_s + _a_n + 1; ++k) {
        auto jj{j - _pml_h_c_start};
        correctPML(_ha(i, j, k), _ha_psi_eb(i, jj, k), _coeff_a_h(jj),
                   _coeff_b_h(jj), _eb(i, j + 1, k), _eb(i, j, k),
                   _c_ha_psi_eb(i, jj, k));
      }
    }
  }

  for (std::size_t i{_b_s}; i < _b_s + _b_n + 1; ++i) {
    for (std::size_t j{_c_h_s}; j < _c_h_s + _c_h_n; ++j) {
      for (std::size_t k{_a_s}; k < _a_s + _a_n; ++k) {
        auto jj{j - _pml_h_c_start};
        correctPML(_hb(i, j, k), _hb_psi_ea(i, jj, k), _coeff_a_h(jj),
                   _coeff_b_h(jj), _ea(i, j + 1, k), _ea(i, j, k),
                   _c_hb_psi_ea(i, jj, k));
      }
    }
  }
}

void PMLCorrectorZ::correctE() {
  for (std::size_t i{_a_s}; i < _a_s + _a_n; ++i) {
    for (std::size_t j{_b_s}; j < _b_s + _b_n + 1; ++j) {
      for (std::size_t k{_c_e_s}; k < _c_e_s + _c_e_n; ++k) {
        auto kk{k - _pml_e_c_start};
        correctPML(_ea(i, j, k), _ea_psi_hb(i, j, kk), _coeff_a_e(kk),
                   _coeff_b_e(kk), _hb(i, j, k), _hb(i, j, k - 1),
                   _c_ea_psi_hb(i, j, kk));
      }
    }
  }
  for (std::size_t i{_a_s}; i < _a_s + _a_n + 1; ++i) {
    for (std::size_t j{_b_s}; j < _b_s + _b_n; ++j) {
      for (std::size_t k{_c_e_s}; k < _c_e_s + _c_e_n; ++k) {
        auto kk{k - _pml_e_c_start};
        correctPML(_eb(i, j, k), _eb_psi_ha(i, j, kk), _coeff_a_e(kk),
                   _coeff_b_e(kk), _ha(i, j, k), _ha(i, j, k - 1),
                   _c_eb_psi_ha(i, j, kk));
      }
    }
  }
}

void PMLCorrectorZ::correctH() {
  for (std::size_t i{_a_s}; i < _a_s + _a_n + 1; ++i) {
    for (std::size_t j{_b_s}; j < _b_s + _b_n; ++j) {
      for (std::size_t k{_c_h_s}; k < _c_h_s + _c_h_n; ++k) {
        auto kk{k - _pml_h_c_start};
        correctPML(_ha(i, j, k), _ha_psi_eb(i, j, kk), _coeff_a_h(kk),
                   _coeff_b_h(kk), _eb(i, j, k + 1), _eb(i, j, k),
                   _c_ha_psi_eb(i, j, kk));
      }
    }
  }
  for (std::size_t i{_a_s}; i < _a_s + _a_n; ++i) {
    for (std::size_t j{_b_s}; j < _b_s + _b_n + 1; ++j) {
      for (std::size_t k{_c_h_s}; k < _c_h_s + _c_h_n; ++k) {
        auto kk{k - _pml_h_c_start};
        correctPML(_hb(i, j, k), _hb_psi_ea(i, j, kk), _coeff_a_h(kk),
                   _coeff_b_h(kk), _ea(i, j, k + 1), _ea(i, j, k),
                   _c_hb_psi_ea(i, j, kk));
      }
    }
  }
}

}  // namespace xfdtd
