#include "waveform_source/tfsf_corrector.h"

#include <sstream>

#include "xfdtd/common/type_define.h"
#include "xfdtd/coordinate_system/coordinate_system.h"
#include "xfdtd/electromagnetic_field/electromagnetic_field.h"

namespace xfdtd {

std::string TFSFCorrector::toString() const {
  std::stringstream ss;
  return ss.str();
}

Real TFSFCorrector::exInc(Index t, Index i, Index j, Index k) const {
  i = i - globalStartI() + 1;
  j = j - globalStartJ();
  k = k - globalStartK();
  const auto& projection_x_half = *_projection_x_half;
  const auto& projection_y_int = *_projection_y_int;
  const auto& projection_z_int = *_projection_z_int;

  auto projection{projection_x_half(i) + projection_y_int(j) +
                  projection_z_int(k)};
  auto index{static_cast<Index>(projection)};
  auto weight{projection - index};
  const auto& e_inc = *_e_inc;

  return _transform_e_x *
         ((1 - weight) * e_inc(t, index) + weight * e_inc(t, index + 1));
}

Real TFSFCorrector::eyInc(Index t, Index i, Index j, Index k) const {
  i = i - globalStartI();
  j = j - globalStartJ() + 1;
  k = k - globalStartK();
  const auto& projection_x_int = *_projection_x_int;
  const auto& projection_y_half = *_projection_y_half;
  const auto& projection_z_int = *_projection_z_int;

  auto projection{projection_x_int(i) + projection_y_half(j) +
                  projection_z_int(k)};
  auto index{static_cast<std::size_t>(projection)};
  auto weight{projection - index};
  const auto& e_inc = *_e_inc;

  return _transform_e_y *
         ((1 - weight) * e_inc(t, index) + weight * e_inc(t, index + 1));
}

Real TFSFCorrector::ezInc(Index t, Index i, Index j, Index k) const {
  i = i - globalStartI();
  j = j - globalStartJ();
  k = k - globalStartK() + 1;
  const auto& projection_x_int = *_projection_x_int;
  const auto& projection_y_int = *_projection_y_int;
  const auto& projection_z_half = *_projection_z_half;

  auto projection{projection_x_int(i) + projection_y_int(j) +
                  projection_z_half(k)};
  auto index{static_cast<std::size_t>(projection)};
  auto weight{projection - index};
  const auto& e_inc = *_e_inc;

  return _transform_e_z *
         ((1 - weight) * e_inc.at(t, index) + weight * e_inc.at(t, index + 1));
}

Real TFSFCorrector::hxInc(Index t, Index i, Index j, Index k) const {
  i = i - globalStartI();
  j = j - globalStartJ() + 1;
  k = k - globalStartK() + 1;
  const auto& projection_x_int = *_projection_x_int;
  const auto& projection_y_half = *_projection_y_half;
  const auto& projection_z_half = *_projection_z_half;

  auto projection{projection_x_int(i) + projection_y_half(j) +
                  projection_z_half(k) - 0.5};
  auto index{static_cast<std::size_t>(projection)};
  auto weight{projection - index};
  const auto& h_inc = *_h_inc;

  return _transform_h_x *
         ((1 - weight) * h_inc(t, index) + weight * h_inc(t, index + 1));
}

Real TFSFCorrector::hyInc(Index t, Index i, Index j, Index k) const {
  i = i - globalStartI() + 1;
  j = j - globalStartJ();
  k = k - globalStartK() + 1;
  const auto& projection_x_half = *_projection_x_half;
  const auto& projection_y_int = *_projection_y_int;
  const auto& projection_z_half = *_projection_z_half;

  auto projection{projection_x_half(i) + projection_y_int(j) +
                  projection_z_half(k) - 0.5};
  auto index{static_cast<std::size_t>(projection)};
  auto weight{projection - index};
  const auto& h_inc = *_h_inc;

  return _transform_h_y *
         ((1 - weight) * h_inc(t, index) + weight * h_inc(t, index + 1));
}

Real TFSFCorrector::hzInc(Index t, Index i, Index j, Index k) const {
  i = i - globalStartI() + 1;
  j = j - globalStartJ() + 1;
  k = k - globalStartK();
  const auto& projection_x_half = *_projection_x_half;
  const auto& projection_y_half = *_projection_y_half;
  const auto& projection_z_int = *_projection_z_int;

  auto projection{projection_x_half(i) + projection_y_half(j) +
                  projection_z_int(k) - 0.5};
  auto index{static_cast<std::size_t>(projection)};
  auto weight{projection - index};
  const auto& h_inc = *_h_inc;

  return _transform_h_z *
         ((1 - weight) * h_inc(t, index) + weight * h_inc(t, index + 1));
}

void TFSF1DCorrector::correctE() {
  auto task = this->task();
  if (_forward) {
    correctTFSF<Axis::Direction::ZN, EMF::Attribute::E, Axis::XYZ::X>(task, 0,
                                                                      0, 0);
    // correctExZN();
    return;
  }

  correctTFSF<Axis::Direction::ZP, EMF::Attribute::E, Axis::XYZ::X>(
      task, _node_offset_i, _node_offset_j, _node_offset_k);
  // correctExZP();
}

void TFSF1DCorrector::correctH() {
  auto task = this->task();
  if (_forward) {
    correctTFSF<Axis::Direction::ZN, EMF::Attribute::H, Axis::XYZ::Y>(
        task, _node_offset_i, _node_offset_j, _node_offset_k);

    // correctHyZN();
    return;
  }

  correctTFSF<Axis::Direction::ZP, EMF::Attribute::H, Axis::XYZ::Y>(
      task, _node_offset_i, _node_offset_j, _node_offset_k);
  // correctHyZP();
}

void TFSF2DCorrector::correctE() {
  auto task = this->task();
  correctTFSF<Axis::Direction::XN, EMF::Attribute::E, Axis::XYZ::Z>(
      task, _node_offset_i, _node_offset_j, _node_offset_k);
  correctTFSF<Axis::Direction::XP, EMF::Attribute::E, Axis::XYZ::Z>(
      task, _node_offset_i, _node_offset_j, _node_offset_k);
  correctTFSF<Axis::Direction::YN, EMF::Attribute::E, Axis::XYZ::Z>(
      task, _node_offset_i, _node_offset_j, _node_offset_k);
  correctTFSF<Axis::Direction::YP, EMF::Attribute::E, Axis::XYZ::Z>(
      task, _node_offset_i, _node_offset_j, _node_offset_k);
}

void TFSF2DCorrector::correctH() {
  auto task = this->task();
  correctTFSF<Axis::Direction::XN, EMF::Attribute::H, Axis::XYZ::Y>(
      task, _node_offset_i, _node_offset_j, _node_offset_k);
  correctTFSF<Axis::Direction::XP, EMF::Attribute::H, Axis::XYZ::Y>(
      task, _node_offset_i, _node_offset_j, _node_offset_k);
  correctTFSF<Axis::Direction::YN, EMF::Attribute::H, Axis::XYZ::X>(
      task, _node_offset_i, _node_offset_j, _node_offset_k);
  correctTFSF<Axis::Direction::YP, EMF::Attribute::H, Axis::XYZ::X>(
      task, _node_offset_i, _node_offset_j, _node_offset_k);
}

void TFSF3DCorrector::correctE() {
  auto task = this->task();
  correctTFSF<Axis::Direction::XN, EMF::Attribute::E, Axis::XYZ::Z>(
      task, _node_offset_i, _node_offset_j, _node_offset_k);
  correctTFSF<Axis::Direction::XP, EMF::Attribute::E, Axis::XYZ::Z>(
      task, _node_offset_i, _node_offset_j, _node_offset_k);
  correctTFSF<Axis::Direction::YN, EMF::Attribute::E, Axis::XYZ::Z>(
      task, _node_offset_i, _node_offset_j, _node_offset_k);
  correctTFSF<Axis::Direction::YP, EMF::Attribute::E, Axis::XYZ::Z>(
      task, _node_offset_i, _node_offset_j, _node_offset_k);
  correctTFSF<Axis::Direction::ZN, EMF::Attribute::E, Axis::XYZ::Y>(
      task, _node_offset_i, _node_offset_j, _node_offset_k);
  correctTFSF<Axis::Direction::ZP, EMF::Attribute::E, Axis::XYZ::Y>(
      task, _node_offset_i, _node_offset_j, _node_offset_k);
  correctTFSF<Axis::Direction::XN, EMF::Attribute::E, Axis::XYZ::Y>(
      task, _node_offset_i, _node_offset_j, _node_offset_k);
  correctTFSF<Axis::Direction::XP, EMF::Attribute::E, Axis::XYZ::Y>(
      task, _node_offset_i, _node_offset_j, _node_offset_k);
  correctTFSF<Axis::Direction::ZN, EMF::Attribute::E, Axis::XYZ::X>(
      task, _node_offset_i, _node_offset_j, _node_offset_k);
  correctTFSF<Axis::Direction::ZP, EMF::Attribute::E, Axis::XYZ::X>(
      task, _node_offset_i, _node_offset_j, _node_offset_k);
  correctTFSF<Axis::Direction::YN, EMF::Attribute::E, Axis::XYZ::X>(
      task, _node_offset_i, _node_offset_j, _node_offset_k);
  correctTFSF<Axis::Direction::YP, EMF::Attribute::E, Axis::XYZ::X>(
      task, _node_offset_i, _node_offset_j, _node_offset_k);
}

void TFSF3DCorrector::correctH() {
  auto task = this->task();

  correctTFSF<Axis::Direction::XN, EMF::Attribute::H, Axis::XYZ::Y>(
      task, _node_offset_i, _node_offset_j, _node_offset_k);
  correctTFSF<Axis::Direction::XP, EMF::Attribute::H, Axis::XYZ::Y>(
      task, _node_offset_i, _node_offset_j, _node_offset_k);
  correctTFSF<Axis::Direction::YN, EMF::Attribute::H, Axis::XYZ::X>(
      task, _node_offset_i, _node_offset_j, _node_offset_k);
  correctTFSF<Axis::Direction::YP, EMF::Attribute::H, Axis::XYZ::X>(
      task, _node_offset_i, _node_offset_j, _node_offset_k);
  correctTFSF<Axis::Direction::ZN, EMF::Attribute::H, Axis::XYZ::X>(
      task, _node_offset_i, _node_offset_j, _node_offset_k);
  correctTFSF<Axis::Direction::ZP, EMF::Attribute::H, Axis::XYZ::X>(
      task, _node_offset_i, _node_offset_j, _node_offset_k);
  correctTFSF<Axis::Direction::XN, EMF::Attribute::H, Axis::XYZ::Z>(
      task, _node_offset_i, _node_offset_j, _node_offset_k);
  correctTFSF<Axis::Direction::XP, EMF::Attribute::H, Axis::XYZ::Z>(
      task, _node_offset_i, _node_offset_j, _node_offset_k);
  correctTFSF<Axis::Direction::YN, EMF::Attribute::H, Axis::XYZ::Z>(
      task, _node_offset_i, _node_offset_j, _node_offset_k);
  correctTFSF<Axis::Direction::YP, EMF::Attribute::H, Axis::XYZ::Z>(
      task, _node_offset_i, _node_offset_j, _node_offset_k);
  correctTFSF<Axis::Direction::ZN, EMF::Attribute::H, Axis::XYZ::Y>(
      task, _node_offset_i, _node_offset_j, _node_offset_k);
  correctTFSF<Axis::Direction::ZP, EMF::Attribute::H, Axis::XYZ::Y>(
      task, _node_offset_i, _node_offset_j, _node_offset_k);
}

}  // namespace xfdtd
