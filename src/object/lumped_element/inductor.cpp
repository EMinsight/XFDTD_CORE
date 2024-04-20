#include <xfdtd/object/lumped_element/inductor.h>

#include <memory>
#include <xtensor.hpp>

#include "corrector/corrector.h"
#include "object/lumped_element_corrector.h"
#include "xfdtd/coordinate_system/coordinate_system.h"

namespace xfdtd {

Inductor::Inductor(std::string name, std::unique_ptr<Cube> cube, Axis::XYZ xyz,
                   Real inductance, std::unique_ptr<Material> material)
    : LumpedElement{std::move(name), std::move(cube), xyz, std::move(material)},
      _inductance{inductance} {
  if (_inductance == 0) {
    _inductance = 1e-10;
  }
}

std::string Inductor::toString() const { return ""; }

void Inductor::init(std::shared_ptr<const GridSpace> grid_space,
                    std::shared_ptr<CalculationParam> calculation_param,
                    std::shared_ptr<EMF> emf) {
  LumpedElement::init(std::move(grid_space), std::move(calculation_param),
                      std::move(emf));

  auto in_f{[](Real i, std::size_t na, std::size_t nb, std::size_t nc) {
    return i * na * nb / nc;
  }};

  auto dx_dy_dz{[](const auto& x, const auto& y, const auto& z, auto&& x_range,
                   auto&& y_range, auto&& z_range) {
    return xt::meshgrid(xt::view(x, x_range), xt::view(y, y_range),
                        xt::view(z, z_range));
  }};

  auto dt{calculationParamPtr()->timeParam()->dt()};
  auto g{gridSpacePtr()};

  _inductance_factor = in_f(_inductance, globalCountSubAxisA(),
                            globalCountSubAxisB(), globalCountMainAxis());

  if (xyz() == Axis::XYZ::X) {
    auto [dx, dy, dz] = dx_dy_dz(g->eSizeX(), g->hSizeY(), g->hSizeZ(),
                                 rangeX(), rangeY(), rangeZ());

    _da = dy;
    _db = dz;
    _dc = dx;
  }

  if (xyz() == Axis::XYZ::Y) {
    auto [dx, dy, dz] = dx_dy_dz(g->hSizeX(), g->eSizeY(), g->hSizeZ(),
                                 rangeX(), rangeY(), rangeZ());

    _da = dz;
    _db = dx;
    _dc = dy;
  }

  if (xyz() == Axis::XYZ::Z) {
    auto [dx, dy, dz] = dx_dy_dz(g->hSizeX(), g->hSizeY(), g->eSizeZ(),
                                 rangeX(), rangeY(), rangeZ());

    _da = dx;
    _db = dy;
    _dc = dz;
  }

  _cjcec = dt * _dc / (_inductance_factor * _da * _db);
  _j = xt::zeros<Real>({nodeCountX(), nodeCountY(), nodeCountZ()});
}

void Inductor::correctUpdateCoefficient() {
  auto calc_param{calculationParamPtr()};
  auto dt{calc_param->timeParam()->dt()};
  if (xyz() == Axis::XYZ::X) {
    const auto eps_view{xt::view(calc_param->materialParam()->epsX(), rangeX(),
                                 rangeY(), rangeZ())};
    const auto sigma_view{xt::view(calc_param->materialParam()->sigmaEX(),
                                   rangeX(), rangeY(), rangeZ())};
    _cecjc = -(2 * dt) / (2 * eps_view + dt * sigma_view);
    return;
  }

  if (xyz() == Axis::XYZ::Y) {
    const auto eps_view{xt::view(calc_param->materialParam()->epsY(), rangeX(),
                                 rangeY(), rangeZ())};
    const auto sigma_view{xt::view(calc_param->materialParam()->sigmaEY(),
                                   rangeX(), rangeY(), rangeZ())};
    _cecjc = -(2 * dt) / (2 * eps_view + dt * sigma_view);
    return;
  }

  if (xyz() == Axis::XYZ::Z) {
    const auto eps_view{xt::view(calc_param->materialParam()->epsZ(), rangeX(),
                                 rangeY(), rangeZ())};
    const auto sigma_view{xt::view(calc_param->materialParam()->sigmaEZ(),
                                   rangeX(), rangeY(), rangeZ())};
    _cecjc = -(2 * dt) / (2 * eps_view + dt * sigma_view);
    return;
  }
}

void Inductor::correctE() {}

void Inductor::correctH() {}

std::unique_ptr<Corrector> Inductor::generateCorrector(
    const Task<std::size_t>& task) {
  if (!taskContainLumpedElement(task)) {
    return nullptr;
  }

  auto domain = makeIndexTask();
  auto intersection = taskIntersection(task, domain);

  if (!intersection.has_value()) {
    return nullptr;
  }

  auto local_task = makeTask(
      makeRange(intersection->xRange().start() - domain.xRange().start(),
                intersection->xRange().end() - domain.xRange().start()),
      makeRange(intersection->yRange().start() - domain.yRange().start(),
                intersection->yRange().end() - domain.yRange().start()),
      makeRange(intersection->zRange().start() - domain.zRange().start(),
                intersection->zRange().end() - domain.zRange().start()));

  return std::make_unique<InductorCorrector>(
      intersection.value(), local_task, calculationParam(),
      fieldMainAxis(EMF::Attribute::E), _j, _cecjc, _cjcec);
}

}  // namespace xfdtd
