#include "updator/updator.h"

#include <utility>

#include "divider/divider.h"

namespace xfdtd {

Updator::Updator(std::shared_ptr<const GridSpace> grid_space,
                 std::shared_ptr<const CalculationParam> calculation_param,
                 std::shared_ptr<EMF> emf, Divider::IndexTask task)
    : _grid_space(std::move(grid_space)),
      _calculation_param(std::move(calculation_param)),
      _emf(std::move(emf)),
      _task{std::move(task)} {}

}  // namespace xfdtd
