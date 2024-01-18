#include <xfdtd/boundary/boundary.h>

namespace xfdtd {

void Boundary::defaultInit(std::shared_ptr<const GridSpace> grid_space,
                           std::shared_ptr<CalculationParam> calculation_param,
                           std::shared_ptr<EMF> emf) {
  _grid_space = std::move(grid_space);
  _calculation_param = std::move(calculation_param);
  _emf = std::move(emf);
}

const GridSpace* Boundary::gridSpacePtr() const { return _grid_space.get(); }

CalculationParam* Boundary::calculationParamPtr() const {
  return _calculation_param.get();
}

EMF* Boundary::emfPtr() const { return _emf.get(); }

}  // namespace xfdtd
