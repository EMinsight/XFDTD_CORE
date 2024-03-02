#include <xfdtd/monitor/voltage_monitor.h>

#include <utility>
#include <xtensor.hpp>

#include "xfdtd/coordinate_system/coordinate_system.h"
#include "xfdtd/monitor/monitor.h"

namespace xfdtd {

VoltageMonitor::VoltageMonitor(std::string name, std::unique_ptr<Cube> cube,
                               Axis::Direction direction,
                               std::string output_dir)
    : Monitor{std::move(cube), std::move(name), std::move(output_dir)},
      _direction{direction} {}

void VoltageMonitor::init(
    std::shared_ptr<const GridSpace> grid_space,
    std::shared_ptr<const CalculationParam> calculation_param,
    std::shared_ptr<const EMF> emf) {
  defaultInit(std::move(grid_space), std::move(calculation_param),
              std::move(emf));

  auto grid_box{gridBoxPtr()};
  _is = grid_box->origin().i();
  _ie = grid_box->end().i();
  _js = grid_box->origin().j();
  _je = grid_box->end().j();
  _ks = grid_box->origin().k();
  _ke = grid_box->end().k();

  if (Axis::formDirectionToXYZ(_direction) == Axis::XYZ::X) {
    _dc = xt::view(gridSpacePtr()->hSizeX(), xt::range(_is, _ie));
    _coff = -_dc / ((_je - _js + 1) * (_ke - _ks + 1));
    if (_direction == Axis::Direction::XN) {
      _coff *= -1.0;
    }
  }

  if (Axis::formDirectionToXYZ(_direction) == Axis::XYZ::Y) {
    _dc = xt::view(gridSpacePtr()->hSizeY(), xt::range(_js, _je));
    _coff = -_dc / ((_ie - _is + 1) * (_ke - _ks + 1));
    if (_direction == Axis::Direction::YN) {
      _coff *= -1.0;
    }
  }

  if (Axis::formDirectionToXYZ(_direction) == Axis::XYZ::Z) {
    _dc = xt::view(gridSpacePtr()->hSizeZ(), xt::range(_ks, _ke));
    _coff = -_dc / ((_ie - _is + 1) * (_je - _js + 1));
    if (_direction == Axis::Direction::ZN) {
      _coff *= -1.0;
    }
  }

  data() = xt::zeros<double>({calculationParamPtr()->timeParam()->size()});
  data() = xt::stack(
      xt::xtuple(calculationParamPtr()->timeParam()->eTime(), data()));
}

void VoltageMonitor::update() {
  auto emf{emfPtr()};
  auto t{calculationParamPtr()->timeParam()->currentTimeStep()};

  if (Axis::formDirectionToXYZ(_direction) == Axis::XYZ::X) {
    for (size_t i{_is}; i < _ie; ++i) {
      for (size_t j{_js}; j < _je + 1; ++j) {
        for (size_t k{_ks}; k < _ke + 1; ++k) {
          data()(1, t) += _coff(i - _is) * emf->ex()(i, j, k);
        }
      }
    }
  }

  if (Axis::formDirectionToXYZ(_direction) == Axis::XYZ::Y) {
    for (size_t i{_is}; i < _ie + 1; ++i) {
      for (size_t j{_js}; j < _je; ++j) {
        for (size_t k{_ks}; k < _ke + 1; ++k) {
          data()(1, t) += _coff(j - _js) * emf->ey()(i, j, k);
        }
      }
    }
  }

  if (Axis::formDirectionToXYZ(_direction) == Axis::XYZ::Z) {
    for (size_t i{_is}; i < _ie + 1; ++i) {
      for (size_t j{_js}; j < _je + 1; ++j) {
        for (size_t k{_ks}; k < _ke; ++k) {
          data()(1, t) += _coff(k - _ks) * emf->ez()(i, j, k);
        }
      }
    }
  }
}

}  // namespace xfdtd
