#include "updator/basic_updator.h"

#include <memory>
#include <utility>
#include <xtensor.hpp>

#include "updator/update_scheme.h"
#include "updator/updator.h"
#include "xfdtd/calculation_param/calculation_param.h"
#include "xfdtd/common/type_define.h"
#include "xfdtd/electromagnetic_field/electromagnetic_field.h"
#include "xfdtd/util/fdtd_basic.h"

namespace xfdtd {

BasicUpdator::BasicUpdator(
    std::shared_ptr<const GridSpace> grid_space,
    std::shared_ptr<const CalculationParam> calculation_param,
    std::shared_ptr<EMF> emf, IndexTask task)
    : Updator(std::move(grid_space), std::move(calculation_param),
              std::move(emf), task) {}

void BasicUpdator::updateH() {
  const auto task = this->task();
  const auto x_range = task.xRange();
  const auto y_range = task.yRange();
  const auto z_range = task.zRange();
  const auto is = basic::GridStructure::hFDTDUpdateXStart(x_range.start());
  const auto ie = basic::GridStructure::hFDTDUpdateXEnd(x_range.end());
  const auto js = basic::GridStructure::hFDTDUpdateYStart(y_range.start());
  const auto je = basic::GridStructure::hFDTDUpdateYEnd(y_range.end());
  const auto ks = basic::GridStructure::hFDTDUpdateZStart(z_range.start());
  const auto ke = basic::GridStructure::hFDTDUpdateZEnd(z_range.end());

  update<EMF::Attribute::H, Axis::XYZ::X>(
      *_emf, *_calculation_param->fdtdCoefficient(), is, ie, js, je, ks, ke);
  update<EMF::Attribute::H, Axis::XYZ::Y>(
      *_emf, *_calculation_param->fdtdCoefficient(), is, ie, js, je, ks, ke);
  update<EMF::Attribute::H, Axis::XYZ::Z>(
      *_emf, *_calculation_param->fdtdCoefficient(), is, ie, js, je, ks, ke);

  // const auto& chxh{_calculation_param->fdtdCoefficient()->chxh()};
  // const auto& chyh{_calculation_param->fdtdCoefficient()->chyh()};
  // const auto& chzh{_calculation_param->fdtdCoefficient()->chzh()};
  // const auto& chxey{_calculation_param->fdtdCoefficient()->chxey()};
  // const auto& chxez{_calculation_param->fdtdCoefficient()->chxez()};
  // const auto& chyex{_calculation_param->fdtdCoefficient()->chyex()};
  // const auto& chyez{_calculation_param->fdtdCoefficient()->chyez()};
  // const auto& chzex{_calculation_param->fdtdCoefficient()->chzex()};
  // const auto& chzey{_calculation_param->fdtdCoefficient()->chzey()};

  // auto& hx{_emf->hx()};
  // auto& hy{_emf->hy()};
  // auto& hz{_emf->hz()};
  // const auto& ex{_emf->ex()};
  // const auto& ey{_emf->ey()};
  // const auto& ez{_emf->ez()};

  // for (std::size_t i{is}; i < ie; ++i) {
  //   for (std::size_t j{js}; j < je; ++j) {
  //     for (std::size_t k{ks}; k < ke; ++k) {
  //       hx(i, j, k) =
  //           hNextDeparted(chxh(i, j, k), hx(i, j, k), chxey(i, j, k), ey(i,
  //           j, k + 1),
  //                 ey(i, j, k), chxez(i, j, k), ez(i, j + 1, k), ez(i, j, k));

  //       hy(i, j, k) =
  //           hNextDeparted(chyh(i, j, k), hy(i, j, k), chyez(i, j, k), ez(i +
  //           1, j, k),
  //                 ez(i, j, k), chyex(i, j, k), ex(i, j, k + 1), ex(i, j, k));

  //       hz(i, j, k) =
  //           hNextDeparted(chzh(i, j, k), hz(i, j, k), chzex(i, j, k), ex(i, j
  //           + 1, k),
  //                 ex(i, j, k), chzey(i, j, k), ey(i + 1, j, k), ey(i, j, k));
  //     }
  //   }
  // }

  // updateHEdge();
}

BasicUpdatorTEM::BasicUpdatorTEM(
    std::shared_ptr<const GridSpace> grid_space,
    std::shared_ptr<const CalculationParam> calculation_param,
    std::shared_ptr<EMF> emf, IndexTask task)
    : BasicUpdator(std::move(grid_space), std::move(calculation_param),
                   std::move(emf), task) {}

void BasicUpdatorTEM::updateE() {
  // const auto ks =
  // basic::GridStructure::exFDTDUpdateZStart(task().zRange().start());
  const auto ks = task().zRange().start() == 0 ? 1 : task().zRange().start();
  const auto ke = basic::GridStructure::exFDTDUpdateZEnd(task().zRange().end());

  const auto& cexe{_calculation_param->fdtdCoefficient()->cexe()};
  const auto& cexhy{_calculation_param->fdtdCoefficient()->cexhy()};

  const auto& hy{_emf->hy()};
  auto& ex{_emf->ex()};

  for (std::size_t k{ks}; k < ke; ++k) {
    ex(0, 0, k) = eNext(cexe(0, 0, k), ex(0, 0, k), cexhy(0, 0, k), hy(0, 0, k),
                        hy(0, 0, k - 1), static_cast<Real>(0.0),
                        static_cast<Real>(0.0), static_cast<Real>(0.0));
  }

  // updateEEdge();
}

// auto BasicUpdatorTEM::updateEEdge() -> void {
//   if (containZNEdge()) {
//     return;
//   }

//   const auto ks = task().zRange().start();

//   const auto& cexe{_calculation_param->fdtdCoefficient()->cexe()};
//   const auto& cexhy{_calculation_param->fdtdCoefficient()->cexhy()};

//   const auto& hy{_emf->hy()};
//   auto& ex{_emf->ex()};
//   auto k = ks;
//   ex(0, 0, k) = eNext(cexe(0, 0, k), ex(0, 0, k), cexhy(0, 0, k), hy(0, 0,
//   k),
//                       hy(0, 0, k - 1), 0.0, 0.0, 0.0);
// }

}  // namespace xfdtd
