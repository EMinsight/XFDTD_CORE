#include "boundary/pml_corrector.h"

#include <xfdtd/common/index_task.h>

#include <sstream>
#include <string>

#include "boundary/pml_scheme.h"

namespace xfdtd {

std::string PMLCorrectorX::toString() const {
  std::stringstream ss;
  ss << "PML X:\n";
  auto e_x_range = makeIndexRange(_c_e_s, _c_e_s + _c_e_n);
  auto h_x_range = makeIndexRange(_c_h_s, _c_h_s + _c_h_n);
  auto y_range = makeIndexRange(_a_s, _a_s + _a_n);
  auto z_range = makeIndexRange(_b_s, _b_s + _b_n);
  ss << " E: " << e_x_range.toString() << ", " << y_range.toString() << ", "
     << z_range.toString() << "\n";
  ss << " H: " << h_x_range.toString() << ", " << y_range.toString() << ", "
     << z_range.toString();
  ss << "\n";
  ss << " PML Global E Start: " << _pml_global_e_start << "\n";
  ss << " PML Global H Start: " << _pml_global_h_start << "\n";
  ss << " Offset C: " << _offset_c << "\n";
  return ss.str();
}

void PMLCorrectorX::correctE() {
  for (std::size_t node_i{_c_e_s}; node_i < _c_e_s + _c_e_n; ++node_i) {
    auto i = node_i - _pml_node_e_start;
    auto global_i = node_i + _offset_c;

    for (std::size_t node_j{_a_s}; node_j < _a_s + _a_n; ++node_j) {
      for (std::size_t node_k{_b_s}; node_k < _b_s + _b_n + 1; ++node_k) {
        correctPML(_ea(node_i, node_j, node_k), _ea_psi_hb(i, node_j, node_k),
                   _coeff_a_e(global_i - _pml_global_e_start),
                   _coeff_b_e(global_i - _pml_global_e_start),
                   _hb(node_i, node_j, node_k), _hb(node_i - 1, node_j, node_k),
                   _c_ea_psi_hb(i, node_j, node_k));
      }
    }

    for (std::size_t node_j{_a_s}; node_j < _a_s + _a_n + 1; ++node_j) {
      for (std::size_t node_k{_b_s}; node_k < _b_s + _b_n; ++node_k) {
        correctPML(_eb(node_i, node_j, node_k), _eb_psi_ha(i, node_j, node_k),
                   _coeff_a_e(global_i - _pml_global_e_start),
                   _coeff_b_e(global_i - _pml_global_e_start),
                   _ha(node_i, node_j, node_k), _ha(node_i - 1, node_j, node_k),
                   _c_eb_psi_ha(i, node_j, node_k));
      }
    }
  }
}

void PMLCorrectorX::correctH() {
  for (std::size_t node_i{_c_h_s}; node_i < _c_h_s + _c_h_n; ++node_i) {
    auto i = node_i - _pml_node_h_start;
    auto global_i = node_i + _offset_c;

    for (std::size_t node_j{_a_s}; node_j < _a_s + _a_n + 1; ++node_j) {
      for (std::size_t node_k{_b_s}; node_k < _b_s + _b_n; ++node_k) {
        correctPML(_ha(node_i, node_j, node_k), _ha_psi_eb(i, node_j, node_k),
                   _coeff_a_h(global_i - _pml_global_h_start),
                   _coeff_b_h(global_i - _pml_global_h_start),
                   _eb(node_i + 1, node_j, node_k), _eb(node_i, node_j, node_k),
                   _c_ha_psi_eb(i, node_j, node_k));
      }
    }

    for (std::size_t node_j{_a_s}; node_j < _a_s + _a_n; ++node_j) {
      for (std::size_t node_k{_b_s}; node_k < _b_s + _b_n + 1; ++node_k) {
        correctPML(_hb(node_i, node_j, node_k), _hb_psi_ea(i, node_j, node_k),
                   _coeff_a_h(global_i - _pml_global_h_start),
                   _coeff_b_h(global_i - _pml_global_h_start),
                   _ea(node_i + 1, node_j, node_k), _ea(node_i, node_j, node_k),
                   _c_hb_psi_ea(i, node_j, node_k));
      }
    }
  }
}

std::string PMLCorrectorY::toString() const {
  std::stringstream ss;
  ss << "PML Y:\n";
  auto e_y_range = makeIndexRange(_c_e_s, _c_e_s + _c_e_n);
  auto h_y_range = makeIndexRange(_c_h_s, _c_h_s + _c_h_n);
  auto x_range = makeIndexRange(_b_s, _b_s + _b_n);
  auto z_range = makeIndexRange(_a_s, _a_s + _a_n);
  ss << " E: " << x_range.toString() << ", " << e_y_range.toString() << ", "
     << z_range.toString() << "\n";
  ss << " H: " << x_range.toString() << ", " << h_y_range.toString() << ", "
     << z_range.toString() << "\n";
  ss << " PML Global E Start: " << _pml_global_e_start << "\n";
  ss << " PML Global H Start: " << _pml_global_h_start << "\n";
  ss << " Offset C: " << _offset_c << "\n";

  return ss.str();
}

void PMLCorrectorY::correctE() {
  for (std::size_t node_i{_b_s}; node_i < _b_s + _b_n + 1; ++node_i) {
    for (std::size_t node_j{_c_e_s}; node_j < _c_e_s + _c_e_n; ++node_j) {
      auto j = node_j - _pml_node_e_start;
      auto global_j = node_j + _offset_c;
      for (std::size_t node_k{_a_s}; node_k < _a_s + _a_n; ++node_k) {
        correctPML(_ea(node_i, node_j, node_k), _ea_psi_hb(node_i, j, node_k),
                   _coeff_a_e(global_j - _pml_global_e_start),
                   _coeff_b_e(global_j - _pml_global_e_start),
                   _hb(node_i, node_j, node_k), _hb(node_i, node_j - 1, node_k),
                   _c_ea_psi_hb(node_i, j, node_k));
      }
    }
  }

  for (std::size_t node_i{_b_s}; node_i < _b_s + _b_n; ++node_i) {
    for (std::size_t node_j{_c_e_s}; node_j < _c_e_s + _c_e_n; ++node_j) {
      auto j = node_j - _pml_node_e_start;
      auto global_j = node_j + _offset_c;
      for (std::size_t node_k{_a_s}; node_k < _a_s + _a_n + 1; ++node_k) {
        correctPML(_eb(node_i, node_j, node_k), _eb_psi_ha(node_i, j, node_k),
                   _coeff_a_e(global_j - _pml_global_e_start),
                   _coeff_b_e(global_j - _pml_global_e_start),
                   _ha(node_i, node_j, node_k), _ha(node_i, node_j - 1, node_k),
                   _c_eb_psi_ha(node_i, j, node_k));
      }
    }
  }
}

void PMLCorrectorY::correctH() {
  for (std::size_t node_i{_b_s}; node_i < _b_s + _b_n; ++node_i) {
    for (std::size_t node_j{_c_h_s}; node_j < _c_h_s + _c_h_n; ++node_j) {
      for (std::size_t node_k{_a_s}; node_k < _a_s + _a_n + 1; ++node_k) {
        auto j = node_j - _pml_node_h_start;
        auto global_j = node_j + _offset_c;
        correctPML(_ha(node_i, node_j, node_k), _ha_psi_eb(node_i, j, node_k),
                   _coeff_a_h(global_j - _pml_global_h_start),
                   _coeff_b_h(global_j - _pml_global_h_start),
                   _eb(node_i, node_j + 1, node_k), _eb(node_i, node_j, node_k),
                   _c_ha_psi_eb(node_i, j, node_k));
      }
    }
  }

  for (std::size_t node_i{_b_s}; node_i < _b_s + _b_n + 1; ++node_i) {
    for (std::size_t node_j{_c_h_s}; node_j < _c_h_s + _c_h_n; ++node_j) {
      for (std::size_t node_k{_a_s}; node_k < _a_s + _a_n; ++node_k) {
        auto j = node_j - _pml_node_h_start;
        auto global_j = node_j + _offset_c;
        correctPML(_hb(node_i, node_j, node_k), _hb_psi_ea(node_i, j, node_k),
                   _coeff_a_h(global_j - _pml_global_h_start),
                   _coeff_b_h(global_j - _pml_global_h_start),
                   _ea(node_i, node_j + 1, node_k), _ea(node_i, node_j, node_k),
                   _c_hb_psi_ea(node_i, j, node_k));
      }
    }
  }
}

std::string PMLCorrectorZ::toString() const {
  std::stringstream ss;
  ss << "PML Z:\n";
  auto e_z_range = makeIndexRange(_c_e_s, _c_e_s + _c_e_n);
  auto h_z_range = makeIndexRange(_c_h_s, _c_h_s + _c_h_n);
  auto x_range = makeIndexRange(_a_s, _a_s + _a_n);
  auto y_range = makeIndexRange(_b_s, _b_s + _b_n);
  ss << " E: " << x_range.toString() << ", " << y_range.toString() << ", "
     << e_z_range.toString() << "\n";
  ss << " H: " << x_range.toString() << ", " << y_range.toString() << ", "
     << h_z_range.toString() << "\n";
  ss << " PML Global E Start: " << _pml_global_e_start << "\n";
  ss << " PML Global H Start: " << _pml_global_h_start << "\n";
  ss << " Offset C: " << _offset_c << "\n";
  return ss.str();
}

void PMLCorrectorZ::correctE() {
  for (std::size_t node_i{_a_s}; node_i < _a_s + _a_n; ++node_i) {
    for (std::size_t node_j{_b_s}; node_j < _b_s + _b_n + 1; ++node_j) {
      for (std::size_t node_k{_c_e_s}; node_k < _c_e_s + _c_e_n; ++node_k) {
        auto k = node_k - _pml_node_e_start;
        auto global_k = node_k + _offset_c;
        correctPML(_ea(node_i, node_j, node_k), _ea_psi_hb(node_i, node_j, k),
                   _coeff_a_e(global_k - _pml_global_e_start),
                   _coeff_b_e(global_k - _pml_global_e_start),
                   _hb(node_i, node_j, node_k), _hb(node_i, node_j, node_k - 1),
                   _c_ea_psi_hb(node_i, node_j, k));
      }
    }
  }
  for (std::size_t node_i{_a_s}; node_i < _a_s + _a_n + 1; ++node_i) {
    for (std::size_t node_j{_b_s}; node_j < _b_s + _b_n; ++node_j) {
      for (std::size_t node_k{_c_e_s}; node_k < _c_e_s + _c_e_n; ++node_k) {
        auto k = node_k - _pml_node_e_start;
        auto global_k = node_k + _offset_c;
        correctPML(_eb(node_i, node_j, node_k), _eb_psi_ha(node_i, node_j, k),
                   _coeff_a_e(global_k - _pml_global_e_start),
                   _coeff_b_e(global_k - _pml_global_e_start),
                   _ha(node_i, node_j, node_k), _ha(node_i, node_j, node_k - 1),
                   _c_eb_psi_ha(node_i, node_j, k));
      }
    }
  }
}

void PMLCorrectorZ::correctH() {
  for (std::size_t node_i{_a_s}; node_i < _a_s + _a_n + 1; ++node_i) {
    for (std::size_t node_j{_b_s}; node_j < _b_s + _b_n; ++node_j) {
      for (std::size_t node_k{_c_h_s}; node_k < _c_h_s + _c_h_n; ++node_k) {
        auto k = node_k - _pml_node_h_start;
        auto global_k = node_k + _offset_c;
        correctPML(_ha(node_i, node_j, node_k), _ha_psi_eb(node_i, node_j, k),
                   _coeff_a_h(global_k - _pml_global_h_start),
                   _coeff_b_h(global_k - _pml_global_h_start),
                   _eb(node_i, node_j, node_k + 1), _eb(node_i, node_j, node_k),
                   _c_ha_psi_eb(node_i, node_j, k));
      }
    }
  }
  for (std::size_t node_i{_a_s}; node_i < _a_s + _a_n; ++node_i) {
    for (std::size_t node_j{_b_s}; node_j < _b_s + _b_n + 1; ++node_j) {
      for (std::size_t node_k{_c_h_s}; node_k < _c_h_s + _c_h_n; ++node_k) {
        auto k = node_k - _pml_node_h_start;
        auto global_k = node_k + _offset_c;
        correctPML(_hb(node_i, node_j, node_k), _hb_psi_ea(node_i, node_j, k),
                   _coeff_a_h(global_k - _pml_global_h_start),
                   _coeff_b_h(global_k - _pml_global_h_start),
                   _ea(node_i, node_j, node_k + 1), _ea(node_i, node_j, node_k),
                   _c_hb_psi_ea(node_i, node_j, k));
      }
    }
  }
}

}  // namespace xfdtd
