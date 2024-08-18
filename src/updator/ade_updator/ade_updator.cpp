#include "updator/ade_updator/ade_updator.h"

namespace xfdtd {

ADEUpdator::ADEUpdator(
    std::shared_ptr<const GridSpace> grid_space,
    std::shared_ptr<const CalculationParam> calculation_param,
    std::shared_ptr<EMF> emf, IndexTask task,
    std::shared_ptr<ADEMethodStorage> ade_method_storage)
    : BasicUpdator{grid_space, calculation_param, emf, task},
      _storage{ade_method_storage} {}

}  // namespace xfdtd
