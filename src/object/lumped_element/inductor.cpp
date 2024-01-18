#include <xfdtd/object/lumped_element/inductor.h>

#include <xtensor.hpp>

#include "xfdtd/coordinate_system/coordinate_system.h"

namespace xfdtd {

Inductor::Inductor(std::string name, std::unique_ptr<Cube> cube, Axis::XYZ xyz,
                   double inductance, std::unique_ptr<Material> material)
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

  auto in_f{[](double i, std::size_t na, std::size_t nb, std::size_t nc) {
    return i * na * nb / nc;
  }};

  auto dx_dy_dz{[](const xt::xarray<double>& x, const xt::xarray<double>& y,
                   const xt::xarray<double>& z, auto&& x_range, auto&& y_range,
                   auto&& z_range) {
    return xt::meshgrid(xt::view(x, x_range), xt::view(y, y_range),
                        xt::view(z, z_range));
  }};

  auto dt{calculationParamPtr()->timeParam()->dt()};
  auto g{gridSpacePtr()};
  auto nx{g->sizeX()};
  auto ny{g->sizeY()};
  auto nz{g->sizeZ()};

  _inductance_factor = in_f(_inductance, nodeCountSubAxisA(),
                            nodeCountSubAxisB(), nodeCountMainAxis());

  if (xyz() == Axis::XYZ::X) {
    auto [dx, dy, dz] = dx_dy_dz(g->eSizeX(), g->hSizeY(), g->hSizeZ(),
                                 rangeX(), rangeY(), rangeZ());

    _da = dy;
    _db = dz;
    _dc = dx;

    emfPtr()->allocateJx(nx, ny + 1, nz + 1);
  }

  if (xyz() == Axis::XYZ::Y) {
    auto [dx, dy, dz] = dx_dy_dz(g->hSizeX(), g->eSizeY(), g->hSizeZ(),
                                 rangeX(), rangeY(), rangeZ());

    _da = dz;
    _db = dx;
    _dc = dy;

    emfPtr()->allocateJy(nx + 1, ny, nz + 1);
  }

  if (xyz() == Axis::XYZ::Z) {
    auto [dx, dy, dz] = dx_dy_dz(g->hSizeX(), g->hSizeY(), g->eSizeZ(),
                                 rangeX(), rangeY(), rangeZ());

    _da = dx;
    _db = dy;
    _dc = dz;

    emfPtr()->allocateJz(nx + 1, ny + 1, nz);
  }

  _cjcec = dt * _dc / (_inductance_factor * _da * _db);
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

void Inductor::correctE() {
  if (xyz() == Axis::XYZ::X) {
    auto& jx{emfPtr()->jx()};
    auto& ex{emfPtr()->ex()};
    auto jx_view{xt::view(jx, rangeX(), rangeY(), rangeZ())};
    auto ex_view{xt::view(ex, rangeX(), rangeY(), rangeZ())};
    ex_view += _cecjc * jx_view;
    jx_view = (_cjcec * ex_view + jx_view);
  }

  if (xyz() == Axis::XYZ::Y) {
    auto& jy{emfPtr()->jy()};
    auto& ey{emfPtr()->ey()};
    auto jy_view{xt::view(jy, rangeX(), rangeY(), rangeZ())};
    auto ey_view{xt::view(ey, rangeX(), rangeY(), rangeZ())};
    ey_view += _cecjc * jy_view;
    jy_view = (_cjcec * ey_view + jy_view);
  }

  if (xyz() == Axis::XYZ::Z) {
    auto& jz{emfPtr()->jz()};
    auto& ez{emfPtr()->ez()};
    auto jz_view{xt::view(jz, rangeX(), rangeY(), rangeZ())};
    auto ez_view{xt::view(ez, rangeX(), rangeY(), rangeZ())};
    ez_view += _cecjc * jz_view;
    jz_view = (_cjcec * ez_view + jz_view);
  }
}

void Inductor::correctH() {}

}  // namespace xfdtd
