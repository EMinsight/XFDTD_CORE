
#include <xfdtd/divider/divider.h>
#include <xfdtd/grid_space/grid_space.h>
#include <xfdtd/monitor/field_monitor.h>
#include <xfdtd/monitor/monitor.h>

#include <utility>
#include <xtensor.hpp>

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

  if (gridSpacePtr()->dimension() == GridSpace::Dimension::ONE) {
    throw XFDTDMonitorException(
        "FieldMonitor cannot be used in 1D simulation(not implemented)");
  }

  auto offset_i = 0;
  auto offset_j = 0;
  auto offset_k = 0;
  if (nodeGridBox().origin().i() == 0) {
    offset_i = 1;
  }
  if (nodeGridBox().origin().j() == 0) {
    offset_j = 1;
  }
  if (nodeGridBox().origin().k() == 0) {
    offset_k =
        gridSpacePtr()->dimension() == GridSpace::Dimension::THREE ? 1 : 0;
  }

  auto new_global_box =
      GridBox{globalGridBox().origin() + Grid{static_cast<size_t>(offset_i),
                                              static_cast<size_t>(offset_j),
                                              static_cast<size_t>(offset_k)},
              globalGridBox().size() - Grid{static_cast<size_t>(offset_i),
                                            static_cast<size_t>(offset_j),
                                            static_cast<size_t>(offset_k)}};

  globalGridBox() = new_global_box;

  auto new_box =
      GridBox{nodeGridBox().origin() + Grid{static_cast<size_t>(offset_i),
                                            static_cast<size_t>(offset_j),
                                            static_cast<size_t>(offset_k)},
              nodeGridBox().size() - Grid{static_cast<size_t>(offset_i),
                                          static_cast<size_t>(offset_j),
                                          static_cast<size_t>(offset_k)}};
  nodeGridBox() = new_box;

  nodeTask() =
      Divider::makeIndexTask(Divider::makeIndexRange(nodeGridBox().origin().i(),
                                                     nodeGridBox().end().i()),
                             Divider::makeIndexRange(nodeGridBox().origin().j(),
                                                     nodeGridBox().end().j()),
                             Divider::makeIndexRange(nodeGridBox().origin().k(),
                                                     nodeGridBox().end().k()));
}

void FieldMonitor::update() {}

void FieldMonitor::output() {
  if (!nodeTask().valid()) {
    return;
  }
  auto em_field{emfPtr()};
  auto calculation_param{calculationParamPtr()};

  auto grid_box{nodeGridBox()};
  auto grid_box_size{grid_box.size()};
  auto x_range{xt::range(grid_box.origin().i(), grid_box.end().i())};
  auto y_range{xt::range(grid_box.origin().j(), grid_box.end().j())};
  auto z_range{xt::range(grid_box.origin().k(), grid_box.end().k())};
  data() = xt::view(em_field->field(field()), x_range, y_range, z_range);
  gatherData();
  Monitor::output();
}

Axis::XYZ FieldMonitor::axis() const { return _axis; }

EMF::Field FieldMonitor::field() const { return _field; }

}  // namespace xfdtd
