#ifndef _XFDTD_CORE_PML_CORRECTOR_H_
#define _XFDTD_CORE_PML_CORRECTOR_H_

#include <xfdtd/common/index_task.h>
#include <xfdtd/common/type_define.h>
#include <xfdtd/coordinate_system/coordinate_system.h>
#include <xfdtd/electromagnetic_field/electromagnetic_field.h>
#include <xfdtd/util/transform/abc_xyz.h>

#include "corrector/corrector.h"

namespace xfdtd {

template <EMF::Attribute attribute, Axis::XYZ xyz>
inline auto correctPML(auto&& field, auto&& psi, const auto& coeff_a,
                       const auto& coeff_b, const auto& dual_field,
                       const auto& c_psi, Index is, Index ie, Index js,
                       Index je, Index ks, Index ke, Index pml_global_start,
                       Index pml_node_start, Index offset_c) {
  for (Index i = is; i < ie; ++i) {
    for (Index j = js; j < je; ++j) {
      for (Index k = ks; k < ke; ++k) {
        auto [a, b, c] = transform::xYZToABC<Index, xyz>(i, j, k);

        const auto global_c = c + offset_c;
        const auto coeff_a_v = coeff_a(global_c - pml_global_start);
        const auto coeff_b_v = coeff_b(global_c - pml_global_start);

        const auto layer_index = c - pml_node_start;
        auto [i_l, j_l, k_l] =
            transform::aBCToXYZ<Index, xyz>(a, b, layer_index);
        const auto c_psi_v = c_psi(i_l, j_l, k_l);
        auto&& psi_v = psi(i_l, j_l, k_l);

        auto&& f_v = field(i, j, k);

        const auto dual_f_p_v = dual_field(i, j, k);
        constexpr auto offset = attribute == EMF::Attribute::E ? -1 : 1;
        auto [i_dual, j_dual, k_dual] =
            transform::aBCToXYZ<Index, xyz>(a, b, c + offset);
        const auto dual_f_q_v = dual_field(i_dual, j_dual, k_dual);

        if constexpr (attribute == EMF::Attribute::E) {
          psi_v = coeff_b_v * psi_v + coeff_a_v * (dual_f_p_v - dual_f_q_v);
        } else {
          psi_v = coeff_b_v * psi_v + coeff_a_v * (dual_f_q_v - dual_f_p_v);
        }
        f_v += c_psi_v * psi_v;
      }
    }
  }
}

template <Axis::XYZ xyz>
class PMLCorrector : public Corrector {
 public:
  PMLCorrector(EMF* emf, IndexTask task, IndexTask node_task,
               Index pml_global_e_start, Index pml_global_h_start,
               Index pml_node_e_start, Index pml_node_h_start, Index offset_c,
               Array1D<Real>* coeff_a_e, Array1D<Real>* coeff_b_e,
               Array1D<Real>* coeff_a_h, Array1D<Real>* coeff_b_h,
               Array3D<Real>* c_ea_psi_hb, Array3D<Real>* c_eb_psi_ha,
               Array3D<Real>* c_ha_psi_eb, Array3D<Real>* c_hb_psi_ea,
               Array3D<Real>* ea_psi_hb, Array3D<Real>* eb_psi_ha,
               Array3D<Real>* ha_psi_eb, Array3D<Real>* hb_psi_ea)
      : _emf{emf},
        _task{task},
        _node_task{node_task},
        _pml_global_e_start{pml_global_e_start},
        _pml_global_h_start{pml_global_h_start},
        _pml_node_e_start{pml_node_e_start},
        _pml_node_h_start{pml_node_h_start},
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
        _hb_psi_ea{hb_psi_ea} {}

  ~PMLCorrector() = default;

  auto correctE() -> void override;

  auto correctH() -> void override;

  auto task() const { return _task; }

  auto nodeTask() const -> IndexTask { return _node_task; }

  auto toString() const -> std::string override {
    std::stringstream ss;
    ss << "PMLCorrector " << Axis::toString(xyz) << " " << task().toString();
    return ss.str();
  }

 private:
  EMF* _emf;
  IndexTask _task, _node_task;

  Index _pml_global_e_start, _pml_global_h_start;
  Index _pml_node_e_start, _pml_node_h_start;
  Index _offset_c;

  Array1D<Real>*_coeff_a_e, *_coeff_b_e, *_coeff_a_h, *_coeff_b_h;
  Array3D<Real>*_c_ea_psi_hb, *_c_eb_psi_ha, *_c_ha_psi_eb, *_c_hb_psi_ea;
  Array3D<Real>*_ea_psi_hb, *_eb_psi_ha, *_ha_psi_eb, *_hb_psi_ea;

  template <EMF::Attribute attribute>
  auto coeffA() const -> Array1D<Real>& {
    if constexpr (attribute == EMF::Attribute::E) {
      return *_coeff_a_e;
    } else {
      return *_coeff_a_h;
    }
  }

  template <EMF::Attribute attribute>
  auto coeffB() const -> Array1D<Real>& {
    if constexpr (attribute == EMF::Attribute::E) {
      return *_coeff_b_e;
    } else {
      return *_coeff_b_h;
    }
  }

  template <EMF::Attribute attribute, Axis::XYZ xyz_0>
  auto cPsi() const -> Array3D<Real>& {
    constexpr auto xyz_a = Axis::tangentialAAxis<xyz>();
    constexpr auto xyz_b = Axis::tangentialBAxis<xyz>();

    static_assert(xyz_0 != xyz, "xyz_0 == xyz!");

    if constexpr (attribute == EMF::Attribute::E) {
      if constexpr (xyz_0 == xyz_a) {
        return *_c_ea_psi_hb;
      } else {
        return *_c_eb_psi_ha;
      }
    } else {
      if constexpr (xyz_0 == xyz_a) {
        return *_c_ha_psi_eb;
      } else {
        return *_c_hb_psi_ea;
      }
    }
  }

  template <EMF::Attribute attribute, Axis::XYZ xyz_0>
  auto psi() -> Array3D<Real>& {
    constexpr auto xyz_a = Axis::tangentialAAxis<xyz>();
    constexpr auto xyz_b = Axis::tangentialBAxis<xyz>();

    static_assert(xyz_0 != xyz, "xyz_0 == xyz!");

    if constexpr (attribute == EMF::Attribute::E) {
      if constexpr (xyz_0 == xyz_a) {
        return *_ea_psi_hb;
      } else {
        return *_eb_psi_ha;
      }
    } else {
      if constexpr (xyz_0 == xyz_a) {
        return *_ha_psi_eb;
      } else {
        return *_hb_psi_ea;
      }
    }
  }
};

template <Axis::XYZ xyz>
inline auto PMLCorrector<xyz>::correctE() -> void {
  constexpr auto attribute = EMF::Attribute::E;
  constexpr auto dual_attribute = EMF::dualAttribute(attribute);
  constexpr auto xyz_a = Axis::tangentialAAxis<xyz>();
  constexpr auto xyz_b = Axis::tangentialBAxis<xyz>();
  const auto pml_global_start = _pml_global_e_start;
  const auto pml_node_start = _pml_node_e_start;
  const auto task = this->task();
  if (!task.valid()) {
    return;
  }

  auto is = task.xRange().start();
  auto js = task.yRange().start();
  auto ks = task.zRange().start();

  auto main_axis_offset = 0;

  {
    // c == 0?
    auto [a, b, c] = transform::xYZToABC<Index, xyz>(
        task.xRange().start(), task.yRange().start(), task.zRange().start());
    if (c == 0) {
      auto [i, j, k] = transform::aBCToXYZ<Index, xyz>(a, b, c + 1);
      main_axis_offset = 1;
      is = i;
      js = j;
      ks = k;
    }
  }

  const auto offset_c = _offset_c;

  {
    auto&& field = _emf->field<attribute, xyz_a>();
    auto&& psi = this->psi<attribute, xyz_a>();
    const auto& coeff_a = coeffA<attribute>();
    const auto& coeff_b = coeffB<attribute>();
    const auto& dual_field = _emf->field<dual_attribute, xyz_b>();
    const auto& c_psi = cPsi<attribute, xyz_a>();
    // correct EA: [a_s, a_e), [b_s, b_e + 1), [c_s, c_e)
    auto [a, b, c] = transform::xYZToABC<Index, xyz>(nodeTask().xRange().end(),
                                                     nodeTask().yRange().end(),
                                                     nodeTask().zRange().end());
    auto [ae, be, ce] = transform::xYZToABC<Index, xyz>(
        task.xRange().end(), task.yRange().end(), task.zRange().end());

    // auto [ie, je, ke] =
    //     transform::aBCToXYZ<Index, xyz>(a, b + 1, c + main_axis_offset);
    Index ie = 0;
    Index je = 0;
    Index ke = 0;
    if (be == b) {
      auto [i, j, k] =
          transform::aBCToXYZ<Index, xyz>(ae, be + 1, ce + main_axis_offset);
      ie = i;
      je = j;
      ke = k;
    } else {
      auto [i, j, k] =
          transform::aBCToXYZ<Index, xyz>(ae, be, ce + main_axis_offset);
      ie = i;
      je = j;
      ke = k;
    }

    correctPML<attribute, xyz>(field, psi, coeff_a, coeff_b, dual_field, c_psi,
                               is, ie, js, je, ks, ke, pml_global_start,
                               pml_node_start, offset_c);
  }

  {
    auto&& field = _emf->field<attribute, xyz_b>();
    auto&& psi = this->psi<attribute, xyz_b>();
    const auto& coeff_a = coeffA<attribute>();
    const auto& coeff_b = coeffB<attribute>();
    const auto& dual_field = _emf->field<dual_attribute, xyz_a>();
    const auto& c_psi = cPsi<attribute, xyz_b>();
    // correct EB: [a_s, a_e + 1), [b_s, b_e), [c_s, c_e)
    auto [a, b, c] = transform::xYZToABC<Index, xyz>(nodeTask().xRange().end(),
                                                     nodeTask().yRange().end(),
                                                     nodeTask().zRange().end());
    auto [ae, be, ce] = transform::xYZToABC<Index, xyz>(
        task.xRange().end(), task.yRange().end(), task.zRange().end());

    // auto [ie, je, ke] =
    //     transform::aBCToXYZ<Index, xyz>(a + 1, b, c + main_axis_offset);

    Index ie = 0;
    Index je = 0;
    Index ke = 0;
    if (ae == a) {
      auto [i, j, k] =
          transform::aBCToXYZ<Index, xyz>(ae + 1, be, ce + main_axis_offset);
      ie = i;
      je = j;
      ke = k;
    } else {
      auto [i, j, k] =
          transform::aBCToXYZ<Index, xyz>(ae, be, ce + main_axis_offset);
      ie = i;
      je = j;
      ke = k;
    }

    correctPML<attribute, xyz>(field, psi, coeff_a, coeff_b, dual_field, c_psi,
                               is, ie, js, je, ks, ke, pml_global_start,
                               pml_node_start, offset_c);
  }
}

template <Axis::XYZ xyz>
inline auto PMLCorrector<xyz>::correctH() -> void {
  constexpr auto attribute = EMF::Attribute::H;
  constexpr auto dual_attribute = EMF::dualAttribute(attribute);
  constexpr auto xyz_a = Axis::tangentialAAxis<xyz>();
  constexpr auto xyz_b = Axis::tangentialBAxis<xyz>();
  const auto pml_global_start = _pml_global_h_start;
  const auto pml_node_start = _pml_node_h_start;
  const auto task = this->task();
  if (!task.valid()) {
    return;
  }

  const auto is = task.xRange().start();
  const auto js = task.yRange().start();
  const auto ks = task.zRange().start();

  const auto offset_c = _offset_c;

  {
    auto&& field = _emf->field<attribute, xyz_a>();
    auto&& psi = this->psi<attribute, xyz_a>();
    const auto& coeff_a = coeffA<attribute>();
    const auto& coeff_b = coeffB<attribute>();
    const auto& dual_field = _emf->field<dual_attribute, xyz_b>();
    const auto& c_psi = cPsi<attribute, xyz_a>();
    // correct HA: [a_s, a_e + 1), [b_s, b_e), [c_s, c_e)
    // auto [a, b, c] = transform::xYZToABC<Index, xyz>(
    //     task().xRange().end(), task().yRange().end(), task().zRange().end());
    // auto [ie, je, ke] = transform::aBCToXYZ<Index, xyz>(a + 1, b, c);
    auto [a, b, c] = transform::xYZToABC<Index, xyz>(nodeTask().xRange().end(),
                                                     nodeTask().yRange().end(),
                                                     nodeTask().zRange().end());
    auto [ae, be, ce] = transform::xYZToABC<Index, xyz>(
        task.xRange().end(), task.yRange().end(), task.zRange().end());

    Index ie = 0;
    Index je = 0;
    Index ke = 0;
    if (ae == a) {
      auto [i, j, k] = transform::aBCToXYZ<Index, xyz>(ae + 1, be, ce);
      ie = i;
      je = j;
      ke = k;
    } else {
      auto [i, j, k] = transform::aBCToXYZ<Index, xyz>(ae, be, ce);
      ie = i;
      je = j;
      ke = k;
    }

    correctPML<attribute, xyz>(field, psi, coeff_a, coeff_b, dual_field, c_psi,
                               is, ie, js, je, ks, ke, pml_global_start,
                               pml_node_start, offset_c);
  }

  {
    auto&& field = _emf->field<attribute, xyz_b>();
    auto&& psi = this->psi<attribute, xyz_b>();
    const auto& coeff_a = coeffA<attribute>();
    const auto& coeff_b = coeffB<attribute>();
    const auto& dual_field = _emf->field<dual_attribute, xyz_a>();
    const auto& c_psi = cPsi<attribute, xyz_b>();
    // correct HB: [a_s, a_e), [b_s, b_e + 1), [c_s, c_e)
    // auto [a, b, c] = transform::xYZToABC<Index, xyz>(
    //     task().xRange().end(), task().yRange().end(), task().zRange().end());
    // auto [ie, je, ke] = transform::aBCToXYZ<Index, xyz>(a, b + 1, c);
    auto [a, b, c] = transform::xYZToABC<Index, xyz>(nodeTask().xRange().end(),
                                                     nodeTask().yRange().end(),
                                                     nodeTask().zRange().end());
    auto [ae, be, ce] = transform::xYZToABC<Index, xyz>(
        task.xRange().end(), task.yRange().end(), task.zRange().end());

    Index ie = 0;
    Index je = 0;
    Index ke = 0;
    if (be == b) {
      auto [i, j, k] = transform::aBCToXYZ<Index, xyz>(ae, be + 1, ce);
      ie = i;
      je = j;
      ke = k;
    } else {
      auto [i, j, k] = transform::aBCToXYZ<Index, xyz>(ae, be, ce);
      ie = i;
      je = j;
      ke = k;
    }

    correctPML<attribute, xyz>(field, psi, coeff_a, coeff_b, dual_field, c_psi,
                               is, ie, js, je, ks, ke, pml_global_start,
                               pml_node_start, offset_c);
  }
}

}  // namespace xfdtd

#endif  // _XFDTD_CORE_PML_CORRECTOR_H_
