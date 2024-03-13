#include "xfdtd/waveform_source/tfsf_1d.h"

#include "xfdtd/util/constant.h"

namespace xfdtd {

TFSF1D::TFSF1D(std::size_t z, bool forward, std::unique_ptr<Waveform> waveform)
    : TFSF{0,
           0,
           z,
           (forward) ? (0) : (constant::PI),
           0,
           0,
           std::move(waveform)},
      _forward{forward} {}

void TFSF1D::init(std::shared_ptr<GridSpace> grid_space,
                  std::shared_ptr<CalculationParam> calculation_param,
                  std::shared_ptr<EMF> emf) {
  TFSF::defaultInit(std::move(grid_space), std::move(calculation_param),
                    std::move(emf));

  // _start = (forward()) ? (boxPtr()->origin().k()) : (boxPtr()->end().k() -
  // 1);
}

std::unique_ptr<Corrector> TFSF1D::generateCorrector(
    const Divider::IndexTask &task) {
  throw std::runtime_error("Not implemented");
}

}  // namespace xfdtd
