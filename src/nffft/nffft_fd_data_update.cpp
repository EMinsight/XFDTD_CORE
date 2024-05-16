#include "nffft/interpolate_scheme.h"
#include "nffft/nffft_fd_data.h"

namespace xfdtd {

auto FDPlaneData::update(std::size_t current_time_step) -> void {
  calculateJ<Axis::Direction::XN>(current_time_step);
  calculateJ<Axis::Direction::XP>(current_time_step);
  calculateJ<Axis::Direction::YN>(current_time_step);
  calculateJ<Axis::Direction::YP>(current_time_step);
  calculateJ<Axis::Direction::ZN>(current_time_step);
  calculateJ<Axis::Direction::ZP>(current_time_step);

  calculateM<Axis::Direction::XN>(current_time_step);
  calculateM<Axis::Direction::XP>(current_time_step);
  calculateM<Axis::Direction::YN>(current_time_step);
  calculateM<Axis::Direction::YP>(current_time_step);
  calculateM<Axis::Direction::ZN>(current_time_step);
  calculateM<Axis::Direction::ZP>(current_time_step);
}

template <Axis::Direction direction>
auto FDPlaneData::calculateM(std::size_t current_time_step) -> void {
  const auto task = this->task<direction>();

  if (!task.valid()) {
    return;
  }

  Real coff_ma = (Axis::directionPositive<direction>()) ? 1.0 : -1.0;
  Real coff_mb = (Axis::directionPositive<direction>()) ? -1.0 : 1.0;

  const auto is = task.xRange().start();
  const auto js = task.yRange().start();
  const auto ks = task.zRange().start();
  const auto ie = task.xRange().end();
  const auto je = task.yRange().end();
  const auto ke = task.zRange().end();

  auto [ma, mb] = surfaceM<direction>();

  constexpr auto xyz = Axis::fromDirectionToXYZ<direction>();
  constexpr auto xyz_a = Axis::tangentialAAxis<xyz>();
  constexpr auto xyz_b = Axis::tangentialBAxis<xyz>();

  constexpr auto attribute = EMF::Attribute::E;
  constexpr auto filed_a =
      EMF::attributeComponentToField(attribute, EMF::xYZToComponent(xyz_a));
  constexpr auto filed_b =
      EMF::attributeComponentToField(attribute, EMF::xYZToComponent(xyz_b));

  const auto emf = _emf.get();

  const auto& ea = emf->field<filed_a>();
  const auto& eb = emf->field<filed_b>();

  for (auto i{is}; i < ie; ++i) {
    for (auto j{js}; j < je; ++j) {
      for (auto k{ks}; k < ke; ++k) {
        ma(i - is, j - js, k - ks) +=
            coff_ma *
            interpolate::interpolateSurfaceCenter<xyz, filed_b>(eb, i, j, k) *
            _transform_e(current_time_step);
        mb(i - is, j - js, k - ks) +=
            coff_mb *
            interpolate::interpolateSurfaceCenter<xyz, filed_a>(ea, i, j, k) *
            _transform_e(current_time_step);
      }
    }
  }
}

template <Axis::Direction direction>
auto FDPlaneData::calculateJ(std::size_t current_time_step) -> void {
  const auto task = this->task<direction>();

  if (!task.valid()) {
    return;
  }

  constexpr Real coff_ja = (Axis::directionPositive<direction>()) ? -1.0 : 1.0;
  constexpr Real coff_jb = (Axis::directionPositive<direction>()) ? 1.0 : -1.0;

  const auto is = task.xRange().start();
  const auto js = task.yRange().start();
  const auto ks = task.zRange().start();
  const auto ie = task.xRange().end();
  const auto je = task.yRange().end();
  const auto ke = task.zRange().end();

  auto [ja, jb] = surfaceJ<direction>();

  constexpr auto xyz = Axis::fromDirectionToXYZ<direction>();
  constexpr auto xyz_a = Axis::tangentialAAxis<xyz>();
  constexpr auto xyz_b = Axis::tangentialBAxis<xyz>();

  constexpr auto attribute = EMF::Attribute::H;
  constexpr auto filed_a =
      EMF::attributeComponentToField(attribute, EMF::xYZToComponent(xyz_a));
  constexpr auto filed_b =
      EMF::attributeComponentToField(attribute, EMF::xYZToComponent(xyz_b));

  const auto emf = _emf.get();

  const auto& ha = emf->field<filed_a>();
  const auto& hb = emf->field<filed_b>();

  for (auto i{is}; i < ie; ++i) {
    for (auto j{js}; j < je; ++j) {
      for (auto k{ks}; k < ke; ++k) {
        ja(i - is, j - js, k - ks) +=
            coff_ja *
            interpolate::interpolateSurfaceCenter<xyz, filed_b>(hb, i, j, k) *
            _transform_h(current_time_step);
        jb(i - is, j - js, k - ks) +=
            coff_jb *
            interpolate::interpolateSurfaceCenter<xyz, filed_a>(ha, i, j, k) *
            _transform_h(current_time_step);
      }
    }
  }
}

}  // namespace xfdtd
