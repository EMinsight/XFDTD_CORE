#ifndef __XFDTD_CORE_NFFFT_TD_DATA_H__
#define __XFDTD_CORE_NFFFT_TD_DATA_H__

#include <xfdtd/calculation_param/calculation_param.h>
#include <xfdtd/common/constant.h>
#include <xfdtd/common/index_task.h>
#include <xfdtd/common/type_define.h>
#include <xfdtd/electromagnetic_field/electromagnetic_field.h>
#include <xfdtd/grid_space/grid_space.h>

#include <algorithm>
#include <memory>
#include <tuple>
#include <utility>
#include <xtensor.hpp>

#include "nffft/interpolate_scheme.h"

namespace xfdtd {

template <Axis::XYZ xyz, EMF::Attribute attribute>
inline auto rVector(Index i, Index j, Index k, const GridSpace* grid_space) {
  if constexpr (xyz == Axis::XYZ::X) {
    return Vector{grid_space->eNodeX()(i), grid_space->hNodeY()(j),
                  grid_space->hNodeZ()(k)};

  } else if constexpr (xyz == Axis::XYZ::Y) {
    return Vector{grid_space->hNodeX()(i), grid_space->eNodeY()(j),
                  grid_space->hNodeZ()(k)};
  } else if constexpr (xyz == Axis::XYZ::Z) {
    return Vector{grid_space->hNodeX()(i), grid_space->hNodeY()(j),
                  grid_space->eNodeZ()(k)};
  }
}

template <Axis::XYZ xyz, EMF::Attribute attribute>
inline auto surfaceArea(Index i, Index j, Index k,
                        const GridSpace* grid_space) -> Real {
  if constexpr (xyz == Axis::XYZ::X) {
    return grid_space->eSizeY()(j) * grid_space->eSizeZ()(k);
  } else if constexpr (xyz == Axis::XYZ::Y) {
    return grid_space->eSizeZ()(k) * grid_space->eSizeX()(i);
  } else if constexpr (xyz == Axis::XYZ::Z) {
    return grid_space->eSizeX()(i) * grid_space->eSizeY()(j);
  }
}

class TDPlaneDataException : public XFDTDException {
 public:
  explicit TDPlaneDataException(const std::string& message)
      : XFDTDException(message) {}
};

template <Axis::Direction direction>
class TDPlaneData {
 public:
  TDPlaneData(std::shared_ptr<const GridSpace> grid_space,
              std::shared_ptr<const CalculationParam> calculation_param,
              std::shared_ptr<const EMF> emf, const IndexTask& task);

  auto task() const -> IndexTask;

  auto initForUpdate(Index num_e, Index num_h, Range<Real> distance_range_e,
                     Range<Real> distance_range_h, Vector r_unit) -> void;

  template <EMF::Attribute attribute>
  auto update() -> void;

  auto wa() const -> Array1D<Real>;

  auto wb() const -> Array1D<Real>;

  auto ua() const -> Array1D<Real>;

  auto ub() const -> Array1D<Real>;

  template <EMF::Attribute attribute>
  auto potential() const
      -> std::tuple<const Array1D<Real>&, const Array1D<Real>&>;

  template <EMF::Attribute attribute>
  auto potential() -> std::tuple<Array1D<Real>&, Array1D<Real>&>;

  template <EMF::Attribute attribute>
  auto distanceRange() const -> Range<Real>;

  std::shared_ptr<const GridSpace> _grid_space;
  std::shared_ptr<const CalculationParam> _calculation_param;
  std::shared_ptr<const EMF> _emf;
  IndexTask _task;
  Vector _r_unit;
  Array2D<Real> _ea_prev, _eb_prev, _ha_prev, _hb_prev;
  Range<Real> _distance_range_e, _distance_range_h;
  Array1D<Real> _wa, _wb, _ua, _ub;

 private:
  auto gridSpacePtr() const -> const GridSpace*;

  auto calculationParamPtr() const -> const CalculationParam*;

  auto emfPrt() const -> const EMF*;

  template <EMF::Attribute attribute>
  auto timeDelay(const Vector& location) const -> Real;

  template <EMF::Attribute attribute>
  auto previousAValue(Index i, Index j, Index k) const -> Real;

  template <EMF::Attribute attribute>
  auto previousBValue(Index i, Index j, Index k) const -> Real;

  template <EMF::Attribute attribute>
  auto setPreviousAValue(Index i, Index j, Index k, Real value) -> void;

  template <EMF::Attribute attribute>
  auto setPreviousBValue(Index i, Index j, Index k, Real value) -> void;
};

template <Axis::XYZ xyz>
inline auto makeArray(const IndexTask& task) {
  if constexpr (xyz == Axis::XYZ::X) {
    return xt::zeros<Real>({task.yRange().size(), task.zRange().size()});
  }

  if constexpr (xyz == Axis::XYZ::Y) {
    return xt::zeros<Real>({task.zRange().size(), task.xRange().size()});
  }

  if constexpr (xyz == Axis::XYZ::Z) {
    return xt::zeros<Real>({task.xRange().size(), task.yRange().size()});
  }

  throw TDPlaneDataException("Invalid Axis::XYZ");
}

template <Axis::Direction D>
TDPlaneData<D>::TDPlaneData(
    std::shared_ptr<const GridSpace> grid_space,
    std::shared_ptr<const CalculationParam> calculation_param,
    std::shared_ptr<const EMF> emf, const IndexTask& task)
    : _grid_space{std::move(grid_space)},
      _calculation_param{std::move(calculation_param)},
      _emf{std::move(emf)},
      _task{task} {
  if (!_task.valid()) {
    return;
  }

  constexpr auto xyz = Axis::fromDirectionToXYZ<D>();

  _ea_prev = makeArray<xyz>(_task);
  _eb_prev = makeArray<xyz>(_task);
  _ha_prev = makeArray<xyz>(_task);
  _hb_prev = makeArray<xyz>(_task);
}

template <Axis::Direction D>
auto TDPlaneData<D>::initForUpdate(Index num_e, Index num_h,
                                   Range<Real> distance_range_e,
                                   Range<Real> distance_range_h,
                                   Vector r_unit) -> void {
  _ua = xt::zeros<Real>({num_e});
  _ub = xt::zeros<Real>({num_e});
  _wa = xt::zeros<Real>({num_h});
  _wb = xt::zeros<Real>({num_h});
  _distance_range_e = distance_range_e;
  _distance_range_h = distance_range_h;
  _r_unit = std::move(r_unit);
}

template <Axis::Direction D>
template <EMF::Attribute attribute>
auto TDPlaneData<D>::update() -> void {
  constexpr auto xyz = Axis::fromDirectionToXYZ<D>();
  constexpr auto xyz_a = Axis::tangentialAAxis<xyz>();
  constexpr auto xyz_b = Axis::tangentialBAxis<xyz>();

  constexpr auto filed_a =
      EMF::attributeComponentToField(attribute, EMF::xYZToComponent(xyz_a));
  constexpr auto filed_b =
      EMF::attributeComponentToField(attribute, EMF::xYZToComponent(xyz_b));

  constexpr Real offset = (attribute == EMF::Attribute::E) ? 0.5 : 0.0;
  Real coeff_a = (attribute == EMF::Attribute::E) ? 1.0 : -1.0;
  if (Axis::directionNegative<D>()) {
    coeff_a *= -1.0;
  }

  Real coeff_b = (attribute == EMF::Attribute::E) ? -1.0 : 1.0;
  if (Axis::directionNegative<D>()) {
    coeff_b *= -1.0;
  }

  auto task = this->task();
  if (!task.valid()) {
    return;
  }

  const EMF* emf = emfPrt();
  auto grid_space = gridSpacePtr();
  auto calculation_param = calculationParamPtr();

  auto current_time_step = calculation_param->timeParam()->currentTimeStep();

  const auto& field_a_value = emf->field<filed_a>();
  const auto& field_b_value = emf->field<filed_b>();

  auto [potential_a, potential_b] = potential<attribute>();

  const auto dt = calculation_param->timeParam()->dt();

  auto is = task.xRange().start();
  auto js = task.yRange().start();
  auto ks = task.zRange().start();
  auto ie = task.xRange().end();
  auto je = task.yRange().end();
  auto ke = task.zRange().end();

  for (auto i{is}; i < ie; ++i) {
    for (auto j{js}; j < je; ++j) {
      for (auto k{ks}; k < ke; ++k) {
        auto r_center = rVector<xyz, attribute>(i, j, k, grid_space);
        auto time_delay = timeDelay<attribute>(r_center);
        auto time_step_delay = time_delay / (constant::C_0 * dt);
        auto nn = static_cast<Index>(
            std::floor(time_step_delay + offset + current_time_step));
        auto coeff_f = time_step_delay + offset + current_time_step - nn;

        auto ds = surfaceArea<xyz, attribute>(i, j, k, grid_space);

        const auto a_avg = interpolate::interpolateSurfaceCenter<xyz, filed_a>(
            field_a_value, i, j, k);
        const auto b_avg = interpolate::interpolateSurfaceCenter<xyz, filed_b>(
            field_b_value, i, j, k);

        auto delta_a =
            (a_avg - previousAValue<attribute>(i - is, j - js, k - ks));
        auto delta_b =
            (b_avg - previousBValue<attribute>(i - is, j - js, k - ks));

        potential_a.at(nn) += coeff_a * (1 - coeff_f) * ds * delta_b;
        potential_a.at(nn + 1) += coeff_a * coeff_f * ds * delta_b;

        potential_b.at(nn) += coeff_b * (1 - coeff_f) * ds * delta_a;
        potential_b.at(nn + 1) += coeff_b * coeff_f * ds * delta_a;

        setPreviousAValue<attribute>(i - is, j - js, k - ks, a_avg);
        setPreviousBValue<attribute>(i - is, j - js, k - ks, b_avg);
      }
    }
  }
}

template <Axis::Direction D>
auto TDPlaneData<D>::task() const -> IndexTask {
  return _task;
}

template <Axis::Direction D>
auto TDPlaneData<D>::wa() const -> Array1D<Real> {
  return _wa / (constant::C_0 * calculationParamPtr()->timeParam()->dt() * 4 *
                constant::PI);
}

template <Axis::Direction D>
auto TDPlaneData<D>::wb() const -> Array1D<Real> {
  return _wb / (constant::C_0 * calculationParamPtr()->timeParam()->dt() * 4 *
                constant::PI);
}

template <Axis::Direction D>
auto TDPlaneData<D>::ua() const -> Array1D<Real> {
  return _ua / (constant::C_0 * calculationParamPtr()->timeParam()->dt() * 4 *
                constant::PI);
}

template <Axis::Direction D>
auto TDPlaneData<D>::ub() const -> Array1D<Real> {
  return _ub / (constant::C_0 * calculationParamPtr()->timeParam()->dt() * 4 *
                constant::PI);
}

template <Axis::Direction D>
template <EMF::Attribute attribute>
auto TDPlaneData<D>::potential() const
    -> std::tuple<const Array1D<Real>&, const Array1D<Real>&> {
  if constexpr (attribute == EMF::Attribute::E) {
    return {_ua, _ub};
  }

  if constexpr (attribute == EMF::Attribute::H) {
    return {_wa, _wb};
  }
}

template <Axis::Direction D>
template <EMF::Attribute attribute>
auto TDPlaneData<D>::potential() -> std::tuple<Array1D<Real>&, Array1D<Real>&> {
  if constexpr (attribute == EMF::Attribute::E) {
    return {_ua, _ub};
  }

  if constexpr (attribute == EMF::Attribute::H) {
    return {_wa, _wb};
  }
}

template <Axis::Direction D>
auto TDPlaneData<D>::gridSpacePtr() const -> const GridSpace* {
  return _grid_space.get();
}

template <Axis::Direction D>
auto TDPlaneData<D>::calculationParamPtr() const -> const CalculationParam* {
  return _calculation_param.get();
}

template <Axis::Direction D>
auto TDPlaneData<D>::emfPrt() const -> const EMF* {
  return _emf.get();
}

template <Axis::Direction D>
template <EMF::Attribute attribute>
auto TDPlaneData<D>::timeDelay(const Vector& location) const -> Real {
  const auto& r_unit = _r_unit;
  if constexpr (attribute == EMF::Attribute::E) {
    return (_distance_range_e.end() -
            (location.x() * r_unit.x() + location.y() * r_unit.y() +
             location.z() * r_unit.z()));
  }

  if constexpr (attribute == EMF::Attribute::H) {
    return (_distance_range_h.end() -
            (location.x() * r_unit.x() + location.y() * r_unit.y() +
             location.z() * r_unit.z()));
  }
}

template <Axis::Direction D>
template <EMF::Attribute attribute>
auto TDPlaneData<D>::distanceRange() const -> Range<Real> {
  auto min_dis = std::numeric_limits<Real>::max();
  auto max_dis = std::numeric_limits<Real>::lowest();

  if (task().valid()) {
    return {min_dis, max_dis};
  }

  const auto& r_unit = _r_unit;

  auto is = task().xRange().start();
  auto js = task().yRange().start();
  auto ks = task().zRange().start();
  auto ie = task().xRange().end();
  auto je = task().yRange().end();
  auto ke = task().zRange().end();

  auto grid_space = gridSpacePtr();

  constexpr auto xyz = Axis::fromDirectionToXYZ<D>();

  for (auto i{is}; i < ie; ++i) {
    for (auto j{js}; j < je; ++j) {
      for (auto k{ks}; k < ke; ++k) {
        auto r = rVector<xyz, attribute>(i, j, k, grid_space);
        auto distance =
            r.x() * r_unit.x() + r.y() * r_unit.y() + r.z() * r_unit.z();
        max_dis = std::max(max_dis, distance);
        min_dis = std::min(min_dis, distance);
      }
    }
  }

  return {min_dis, max_dis};
}

template <Axis::Direction D>
template <EMF::Attribute attribute>
auto TDPlaneData<D>::previousAValue(Index i, Index j, Index k) const -> Real {
  constexpr auto xyz = Axis::fromDirectionToXYZ<D>();
  if constexpr (attribute == EMF::Attribute::E) {
    if constexpr (xyz == Axis::XYZ::X) {
      return _ea_prev.at(j, k);
    } else if constexpr (xyz == Axis::XYZ::Y) {
      return _ea_prev.at(k, i);
    } else if constexpr (xyz == Axis::XYZ::Z) {
      return _ea_prev.at(i, j);
    }
  }

  if constexpr (attribute == EMF::Attribute::H) {
    if constexpr (xyz == Axis::XYZ::X) {
      return _ha_prev.at(j, k);
    } else if constexpr (xyz == Axis::XYZ::Y) {
      return _ha_prev.at(k, i);
    } else if constexpr (xyz == Axis::XYZ::Z) {
      return _ha_prev.at(i, j);
    }
  }
}

template <Axis::Direction D>
template <EMF::Attribute attribute>
auto TDPlaneData<D>::previousBValue(Index i, Index j, Index k) const -> Real {
  constexpr auto xyz = Axis::fromDirectionToXYZ<D>();
  if constexpr (attribute == EMF::Attribute::E) {
    if constexpr (xyz == Axis::XYZ::X) {
      return _eb_prev.at(j, k);
    } else if constexpr (xyz == Axis::XYZ::Y) {
      return _eb_prev.at(k, i);
    } else if constexpr (xyz == Axis::XYZ::Z) {
      return _eb_prev.at(i, j);
    }
  }

  if constexpr (attribute == EMF::Attribute::H) {
    if constexpr (xyz == Axis::XYZ::X) {
      return _hb_prev.at(j, k);
    } else if constexpr (xyz == Axis::XYZ::Y) {
      return _hb_prev.at(k, i);
    } else if constexpr (xyz == Axis::XYZ::Z) {
      return _hb_prev.at(i, j);
    }
  }
}

template <Axis::Direction D>
template <EMF::Attribute attribute>
auto TDPlaneData<D>::setPreviousAValue(Index i, Index j, Index k,
                                       Real value) -> void {
  constexpr auto xyz = Axis::fromDirectionToXYZ<D>();
  if constexpr (attribute == EMF::Attribute::E) {
    if constexpr (xyz == Axis::XYZ::X) {
      _ea_prev.at(j, k) = value;
    } else if constexpr (xyz == Axis::XYZ::Y) {
      _ea_prev.at(k, i) = value;
    } else if constexpr (xyz == Axis::XYZ::Z) {
      _ea_prev.at(i, j) = value;
    }
  }

  if constexpr (attribute == EMF::Attribute::H) {
    if constexpr (xyz == Axis::XYZ::X) {
      _ha_prev.at(j, k) = value;
    } else if constexpr (xyz == Axis::XYZ::Y) {
      _ha_prev.at(k, i) = value;
    } else if constexpr (xyz == Axis::XYZ::Z) {
      _ha_prev.at(i, j) = value;
    }
  }
}

template <Axis::Direction D>
template <EMF::Attribute attribute>
auto TDPlaneData<D>::setPreviousBValue(Index i, Index j, Index k,
                                       Real value) -> void {
  constexpr auto xyz = Axis::fromDirectionToXYZ<D>();
  if constexpr (attribute == EMF::Attribute::E) {
    if constexpr (xyz == Axis::XYZ::X) {
      _eb_prev.at(j, k) = value;
    } else if constexpr (xyz == Axis::XYZ::Y) {
      _eb_prev.at(k, i) = value;
    } else if constexpr (xyz == Axis::XYZ::Z) {
      _eb_prev.at(i, j) = value;
    }
  }

  if constexpr (attribute == EMF::Attribute::H) {
    if constexpr (xyz == Axis::XYZ::X) {
      _hb_prev.at(j, k) = value;
    } else if constexpr (xyz == Axis::XYZ::Y) {
      _hb_prev.at(k, i) = value;
    } else if constexpr (xyz == Axis::XYZ::Z) {
      _hb_prev.at(i, j) = value;
    }
  }
}

}  // namespace xfdtd

#endif  // __XFDTD_CORE_NFFFT_TD_DATA_H__
