#include "updator/updator.h"

#include <sstream>
#include <utility>

#include "divider/divider.h"

namespace xfdtd {

Updator::Updator(std::shared_ptr<const GridSpace> grid_space,
                 std::shared_ptr<const CalculationParam> calculation_param,
                 std::shared_ptr<EMF> emf, Divider::IndexTask task)
    : _grid_space(std::move(grid_space)),
      _calculation_param(std::move(calculation_param)),
      _emf(std::move(emf)),
      _task{task} {}

std::string Updator::toString() const {
  std::stringstream ss;
  ss << "Updator: ";
  ss << _task.toString() << "\n";
  ss << "H field: " << _task.toString();
  return ss.str();
}

}  // namespace xfdtd
