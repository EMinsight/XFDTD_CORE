#include <xfdtd/object/lumped_element/resistor.h>

#include <memory>
#include <utility>
#include <xtensor.hpp>

#include "xfdtd/coordinate_system/coordinate_system.h"

namespace xfdtd {

Resistor::Resistor(std::string name, std::unique_ptr<Cube> cube, Axis::XYZ xyz,
                   double resistance, std::unique_ptr<Material> material)
    : LumpedElement{std::move(name), std::move(cube), xyz, std::move(material)},
      _resistance{resistance} {
  if (_resistance == 0) {
    _resistance = 1e-10;
  }
}

std::string Resistor::toString() const {
  // return "Resistor{name: " + name() + ", shape: " + shape()->toString() +
  //        ", material: " + material()->toString() +
  //        ", direction: " + Axis::toString(_direction) +
  //        ", resistance: " + std::to_string(_resistance) + "}";
  return "";
}

double Resistor::resistance() const { return _resistance; }

void Resistor::init(std::shared_ptr<const GridSpace> grid_space,
                    std::shared_ptr<CalculationParam> calculation_param,
                    std::shared_ptr<EMF> emf) {
  xfdtd::LumpedElement::init(std::move(grid_space),
                             std::move(calculation_param), std::move(emf));

  auto rf{[](double r, std::size_t na, std::size_t nb, std::size_t nc) {
    return r * na * nb / nc;
  }};

  auto dx_dy_dz{[](const xt::xarray<double>& x, const xt::xarray<double>& y,
                   const xt::xarray<double>& z, auto&& x_range, auto&& y_range,
                   auto&& z_range) {
    return xt::meshgrid(xt::view(x, x_range), xt::view(y, y_range),
                        xt::view(z, z_range));
  }};

  _resistance_factor = rf(_resistance, nodeCountSubAxisA(), nodeCountSubAxisB(),
                          nodeCountMainAxis());

  auto g{gridSpacePtr()};

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
}

void Resistor::correctUpdateCoefficient() {
  auto calc_param{calculationParamPtr()};
  auto dt{calc_param->timeParam()->dt()};
  _beta = dt * _dc / (_resistance_factor * _da * _db);
  auto func{[this, dt](auto& cece, auto& cecha, auto& cechb, const auto& eps,
                       const auto& sigma) {
    auto range_x{rangeX()};
    auto range_y{rangeY()};
    auto range_z{rangeZ()};
    auto cece_view{xt::view(cece, range_x, range_y, range_z)};
    auto cecha_view{xt::view(cecha, range_x, range_y, range_z)};
    auto cechb_view{xt::view(cechb, range_x, range_y, range_z)};
    const auto eps_view{xt::view(eps, range_x, range_y, range_z)};
    const auto sigma_view{xt::view(sigma, range_x, range_y, range_z)};

    cece_view = (2 * eps_view - dt * sigma_view - _beta) /
                (2 * eps_view + dt * sigma_view + _beta);
    cecha_view = -2 * dt / ((2 * eps_view + dt * sigma_view + _beta) * _db);
    cechb_view = 2 * dt / ((2 * eps_view + dt * sigma_view + _beta) * _da);
  }};

  if (xyz() == Axis::XYZ::X) {
    func(calc_param->fdtdCoefficient()->cexe(),
         calc_param->fdtdCoefficient()->cexhy(),
         calc_param->fdtdCoefficient()->cexhz(),
         calc_param->materialParam()->epsX(),
         calc_param->materialParam()->sigmaEX());
    return;
  }

  if (xyz() == Axis::XYZ::Y) {
    func(calc_param->fdtdCoefficient()->ceye(),
         calc_param->fdtdCoefficient()->ceyhz(),
         calc_param->fdtdCoefficient()->ceyhx(),
         calc_param->materialParam()->epsY(),
         calc_param->materialParam()->sigmaEY());
    return;
  }

  if (xyz() == Axis::XYZ::Z) {
    func(calc_param->fdtdCoefficient()->ceze(),
         calc_param->fdtdCoefficient()->cezhx(),
         calc_param->fdtdCoefficient()->cezhy(),
         calc_param->materialParam()->epsZ(),
         calc_param->materialParam()->sigmaEZ());
    return;
  }
}

void Resistor::correctE() {}

void Resistor::correctH() {}

}  // namespace xfdtd
