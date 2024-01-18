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

  _start = (forward()) ? (boxPtr()->origin().k()) : (boxPtr()->end().k() - 1);
}

void TFSF1D::correctE() {
  const auto li{boxPtr()->origin().i()};
  const auto lj{boxPtr()->origin().j()};
  const auto lk{boxPtr()->origin().k()};
  const auto ri{boxPtr()->end().i()};
  const auto rj{boxPtr()->end().j()};
  const auto rk{boxPtr()->end().k()};

  // for (std::size_t i{lj}; i < rj; ++i) {
  //   for (std::size_t j{lj}; j < rj + 1; ++j) {
  //     auto &ex_zn{emfPtr()->ex()(i, j, lk)};
  //     auto hy_i{hyInc(i, j, lk - 1)};
  //     ex_zn += hy_i;
  //   }
  // }

  auto &ex_zn{emfPtr()->ex()(0, 0, lk)};
  auto hy_i{hyInc(0, 0, lk - 1)};
  ex_zn += cax() * hy_i;

  auto &ex_zp{emfPtr()->ex()(0, 0, rk)};
  hy_i = {hyInc(0, 0, rk)};
  ex_zp -= cax() * hy_i;
}

void TFSF1D::correctH() {
  const auto li{boxPtr()->origin().i()};
  const auto lj{boxPtr()->origin().j()};
  const auto lk{boxPtr()->origin().k()};
  const auto ri{boxPtr()->end().i()};
  const auto rj{boxPtr()->end().j()};
  const auto rk{boxPtr()->end().k()};

  auto &hy_zn{emfPtr()->hy()(0, 0, lk-1)};
  auto ex_i{exInc(0, 0, lk)};
  hy_zn += cbx() * ex_i;

  auto &hy_zp{emfPtr()->hy()(0, 0, rk)};
  ex_i = exInc(0, 0, rk);
  hy_zp -= cbx() * ex_i;
}

}  // namespace xfdtd
