#include <xfdtd/electromagnetic_field/electromagnetic_field.h>
#include <xfdtd/monitor/field_time_monitor.h>

#include <utility>

#include "xfdtd/grid_space/grid_space.h"

namespace xfdtd {

FieldTimeMonitor::FieldTimeMonitor(std::string name, std::unique_ptr<Cube> cube,
                                   std::string output_dir, EMF::Field field)
    : TimeMonitor(std::move(name), std::move(cube), std::move(output_dir)),
      _field(field) {}

auto FieldTimeMonitor::init(
    std::shared_ptr<const GridSpace> grid_space,
    std::shared_ptr<const CalculationParam> calculation_param,
    std::shared_ptr<const EMF> emf) -> void {
  defaultInit(std::move(grid_space), std::move(calculation_param),
              std::move(emf));

  auto component = EMF::componentFromField(field());
  auto component_to_axis = [](const EMF::Component &c) {
    switch (c) {
      case EMF::Component::X:
        return Axis::XYZ::X;
      case EMF::Component::Y:
        return Axis::XYZ::Y;
      case EMF::Component::Z:
        return Axis::XYZ::Z;
      default:
        throw XFDTDException("Invalid component type");
    }
  };

  auto axis = component_to_axis(component);

  // TODO(franzero): temporary way. need to be refactored
  // Example: The box for Hx is. The box for Ex is
  auto offset_i = 0;
  auto offset_j = 0;
  auto offset_k = 0;
  if (nodeGridBox().origin().i() == 0) {
    offset_i = gridSpacePtr()->dimension() == GridSpace::Dimension::ONE ? 0 : 1;
  }
  if (nodeGridBox().origin().j() == 0) {
    offset_j = gridSpacePtr()->dimension() == GridSpace::Dimension::ONE ? 0 : 1;
  }
  if (nodeGridBox().origin().k() == 0) {
    offset_k =
        gridSpacePtr()->dimension() == GridSpace::Dimension::THREE ? 1 : 0;
  }

  setGlobalGridBox(
      GridBox{globalGridBox().origin() + Grid{static_cast<size_t>(offset_i),
                                              static_cast<size_t>(offset_j),
                                              static_cast<size_t>(offset_k)},
              globalGridBox().size() - Grid{static_cast<size_t>(offset_i),
                                            static_cast<size_t>(offset_j),
                                            static_cast<size_t>(offset_k)}});

  setNodeGridBox(
      GridBox{nodeGridBox().origin() + Grid{static_cast<size_t>(offset_i),
                                            static_cast<size_t>(offset_j),
                                            static_cast<size_t>(offset_k)},
              nodeGridBox().size() - Grid{static_cast<size_t>(offset_i),
                                          static_cast<size_t>(offset_j),
                                          static_cast<size_t>(offset_k)}});

  setGlobalTask(makeIndexTask(
      makeIndexRange(globalGridBox().origin().i(), globalGridBox().end().i()),
      makeIndexRange(globalGridBox().origin().j(), globalGridBox().end().j()),
      makeIndexRange(globalGridBox().origin().k(), globalGridBox().end().k())));

  setNodeTask(makeIndexTask(
      makeIndexRange(nodeGridBox().origin().i(), nodeGridBox().end().i()),
      makeIndexRange(nodeGridBox().origin().j(), nodeGridBox().end().j()),
      makeIndexRange(nodeGridBox().origin().k(), nodeGridBox().end().k())));
}

auto FieldTimeMonitor::update() -> void {
  auto t = calculationParamPtr()->timeParam()->currentTimeStep();
  auto emf = emfPtr();

  for (auto i{nodeTask().xRange().start()}; i < nodeTask().xRange().end();
       ++i) {
    for (auto j{nodeTask().yRange().start()}; j < nodeTask().yRange().end();
         ++j) {
      for (auto k{nodeTask().zRange().start()}; k < nodeTask().zRange().end();
           ++k) {
        data()(t) = emf->field(field())(i, j, k);
      }
    }
  }
}

auto FieldTimeMonitor::initTimeDependentVariable() -> void {
  auto attribute = EMF::attributeFromField(field());
  switch (attribute) {
    case EMF::Attribute::E:
      setTime(calculationParamPtr()->timeParam()->eTime());
      break;
    case EMF::Attribute::H:
      setTime(calculationParamPtr()->timeParam()->hTime());
      break;
    default:
      throw XFDTDEMFException("Invalid attribute type");
  }

  data() = xt::zeros<Real>({time().size()});
}

auto FieldTimeMonitor::initParallelizedConfig() -> void {
  makeMpiSubComm();
  if (valid() && 1 < monitorMpiConfig().size()) {
    throw XFDTDMonitorException(
        "XFDTD Field Time Monitor does not support parallelization");
  }
}

auto FieldTimeMonitor::toString() const -> std::string {
  return "FieldTimeMonitor: " + name();
}

auto FieldTimeMonitor::field() const -> EMF::Field { return _field; }

}  // namespace xfdtd
