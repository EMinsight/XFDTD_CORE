#include "updator/dispersive_material_updator.h"

#include <xfdtd/util/fdtd_basic.h>

#include <memory>
#include <utility>

namespace xfdtd {

LinearDispersiveMaterialADEUpdator::LinearDispersiveMaterialADEUpdator(
    std::shared_ptr<const GridSpace> grid_space,
    std::shared_ptr<const CalculationParam> calculation_param,
    std::shared_ptr<EMF> emf, IndexTask task)
    : BasicUpdator3D(std::move(grid_space), std::move(calculation_param),
                     std::move(emf), task) {
  if (_grid_space->type() != GridSpace::Type::UNIFORM) {
    throw std::runtime_error("Grid space type must be uniform");
  }
}

}  // namespace xfdtd
