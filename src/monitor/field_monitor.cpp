
#include "xfdtd/monitor/field_monitor.h"

#include <utility>
#include <xtensor.hpp>

#include "xfdtd/monitor/monitor.h"

namespace xfdtd {

FieldMonitor::FieldMonitor(std::unique_ptr<Shape> shape, Axis::XYZ axis,
                           EMF::Field field, std::string name,
                           std::string output_dir_path)
    : Monitor{std::move(shape), std::move(name), std::move(output_dir_path)},
      _axis{axis},
      _field{field} {}

void FieldMonitor::init(
    std::shared_ptr<const GridSpace> grid_space,
    std::shared_ptr<const CalculationParam> calculation_param,
    std::shared_ptr<const EMF> emf) {
  defaultInit(grid_space, calculation_param, emf);
}

void FieldMonitor::update() {}

void FieldMonitor::output() {
  auto em_field{emfPtr()};
  auto grid_space{gridSpacePtr()};
  auto calculation_param{calculationParamPtr()};

  auto grid_box{gridBoxPtr()};
  auto grid_box_size{grid_box->size()};
  auto x_range{xt::range(grid_box->origin().i(), grid_box->end().i())};
  auto y_range{xt::range(grid_box->origin().j(), grid_box->end().j())};
  auto z_range{xt::range(grid_box->origin().k(), grid_box->end().k())};
  data() = xt::view(em_field->field(field()), x_range, y_range, z_range);
  Monitor::output();
}

Axis::XYZ FieldMonitor::axis() const { return _axis; }

EMF::Field FieldMonitor::field() const { return _field; }

}  // namespace xfdtd
