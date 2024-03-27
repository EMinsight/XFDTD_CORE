#include <xfdtd/monitor/current_monitor.h>

#include <memory>
#include <xtensor.hpp>

#include "xfdtd/coordinate_system/coordinate_system.h"

namespace xfdtd {

CurrentMonitor::CurrentMonitor(std::string name, std::unique_ptr<Cube> cube,
                               Axis::Direction direction,
                               std::string output_dir)
    : Monitor{std::move(cube), std::move(name), std::move(output_dir)},
      _direction{direction} {}

void CurrentMonitor::init(
    std::shared_ptr<const GridSpace> grid_space,
    std::shared_ptr<const CalculationParam> calculation_param,
    std::shared_ptr<const EMF> emf) {
  defaultInit(grid_space, calculation_param, emf);

  auto grid_box{nodeGridBox()};
  _is = grid_box.origin().i();
  _ie = grid_box.end().i();
  _js = grid_box.origin().j();
  _je = grid_box.end().j();
  _ks = grid_box.origin().k();
  _ke = grid_box.end().k();

  if (Axis::formDirectionToXYZ(_direction) == Axis::XYZ::X) {
    _da = xt::view(gridSpacePtr()->eSizeY(), xt::range(_js, _je + 1));
    _db = xt::view(gridSpacePtr()->eSizeZ(), xt::range(_ks, _ke + 1));
    if (_direction == Axis::Direction::XP) {
      _positive = 1.0;
    } else {
      _positive = -1.0;
    }
  }

  if (Axis::formDirectionToXYZ(_direction) == Axis::XYZ::Y) {
    _da = xt::view(gridSpacePtr()->eSizeZ(), xt::range(_ks, _ke + 1));
    _db = xt::view(gridSpacePtr()->eSizeX(), xt::range(_is, _ie + 1));
    if (_direction == Axis::Direction::YP) {
      _positive = 1.0;
    } else {
      _positive = -1.0;
    }
  }

  if (Axis::formDirectionToXYZ(_direction) == Axis::XYZ::Z) {
    _da = xt::view(gridSpacePtr()->eSizeX(), xt::range(_is, _ie + 1));
    _db = xt::view(gridSpacePtr()->eSizeY(), xt::range(_js, _je + 1));
    if (_direction == Axis::Direction::ZP) {
      _positive = 1.0;
    } else {
      _positive = -1.0;
    }
  }
}

void CurrentMonitor::update() {
  auto emf{emfPtr()};
  auto t{calculationParamPtr()->timeParam()->currentTimeStep()};
  double current{0};
  double integral_x{0};
  double integral_y{0};
  double integral_z{0};

  if (Axis::formDirectionToXYZ(_direction) == Axis::XYZ::X) {
    for (size_t j{_js}; j <= _je; ++j) {
      integral_y +=
          (emf->hy()(_ie - 1, j, _ks - 1) - emf->hy()(_ie - 1, j, _ke)) *
          _da(j - _js);
    }
    for (size_t k{_ks}; k <= _ke; ++k) {
      integral_z +=
          (emf->hz()(_ie - 1, _je, k) - emf->hz()(_ie - 1, _js - 1, k)) *
          _db(k - _ks);
    }
    data()(1, t) = _positive * (integral_z + integral_y);
    return;
  }

  if (Axis::formDirectionToXYZ(_direction) == Axis::XYZ::Y) {
    for (size_t k{_ks}; k <= _ke; ++k) {
      integral_z +=
          (emf->hz()(_is - 1, _je - 1, k) - emf->hz()(_ie, _je - 1, k)) *
          _da(k - _ks);
    }
    for (size_t i{_is}; i <= _ie; ++i) {
      integral_x +=
          (emf->hx()(i, _je - 1, _ke) - emf->hx()(i, _je - 1, _ks - 1)) *
          _db(i - _is);
    }
    data()(1, t) = _positive * (integral_x + integral_z);
    return;
  }

  if (Axis::formDirectionToXYZ(_direction) == Axis::XYZ::Z) {
    for (size_t i{_is}; i <= _ie; ++i) {
      integral_x +=
          (emf->hx()(i, _js - 1, _ke - 1) - emf->hx()(i, _je, _ke - 1)) *
          _da(i - _is);
    }
    for (size_t j{_js}; j <= _je; ++j) {
      integral_y +=
          (emf->hy()(_ie, j, _ke - 1) - emf->hy()(_is - 1, j, _ke - 1)) *
          _db(j - _js);
    }
    data()(1, t) = _positive * (integral_x + integral_y);
    return;
  }
}

void CurrentMonitor::initTimeDependentVariable() {
  data() = xt::zeros<double>({calculationParamPtr()->timeParam()->size()});
  data() = xt::stack(
      xt::xtuple(calculationParamPtr()->timeParam()->hTime(), data()));
}

}  // namespace xfdtd
