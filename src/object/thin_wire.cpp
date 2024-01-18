#include "xfdtd/object/thin_wire.h"

#include <xtensor.hpp>

#include "xfdtd/material/material.h"

namespace xfdtd {

ThinWire::ThinWire(const std::string& name, std::unique_ptr<Cylinder> shape)
    : Object(name, std::move(shape), Material::createAir()) {}

void ThinWire::init(std::shared_ptr<const GridSpace> grid_space,
                    std::shared_ptr<CalculationParam> calculation_param,
                    std::shared_ptr<EMF> emf) {
  Object::init(grid_space, calculation_param, emf);

  auto shape{dynamic_cast<Cylinder*>(shapePtr())};
  if (shape == nullptr) {
    throw XFDTDObjectException("ThinWire::init: shape is not Cylinder");
  }

  if (shape->axis() == Axis::XYZ::X) {
    _axis = Axis::XYZ::X;
  } else if (shape->axis() == Axis::XYZ::Y) {
    _axis = Axis::XYZ::Y;
  } else if (shape->axis() == Axis::XYZ::Z) {
    _axis = Axis::XYZ::Z;
  } else {
    throw XFDTDObjectException("ThinWire::init: invalid axis");
  }
  _radius = shape->radius();

  if (!isEnoughThin()) {
    throw XFDTDObjectException("ThinWire::init: wire is not thin enough");
  }

  if (gridSpacePtr()->type() != GridSpace::Type::UNIFORM) {
    throw XFDTDObjectException("ThinWire::init: grid space is not uniform");
  }
}

void ThinWire::correctUpdateCoefficient() {
  auto box{gridBoxPtr()};
  auto is{box->origin().i()};
  auto ie{box->end().i()};
  auto js{box->origin().j()};
  auto je{box->end().j()};
  auto ks{box->origin().k()};
  auto ke{box->end().k()};

  auto dx{gridSpacePtr()->basedDx()};
  auto dy{gridSpacePtr()->basedDy()};
  auto dz{gridSpacePtr()->basedDz()};

  auto x_range{xt::range(is, ie)};
  auto y_range{xt::range(js, je)};
  auto z_range{xt::range(ks, ke)};

  auto dt{calculationParamPtr()->timeParam()->dt()};

  if (_axis == Axis::XYZ::X) {
    auto cexe_view{xt::view(calculationParamPtr()->fdtdCoefficient()->cexe(),
                            x_range, js, ks)};
    auto cexhy_view{xt::view(calculationParamPtr()->fdtdCoefficient()->cexhy(),
                             x_range, js, ks)};
    auto cexhz_view{xt::view(calculationParamPtr()->fdtdCoefficient()->cexhz(),
                             x_range, js, ks)};
    cexe_view = 0;
    cexhy_view = 0;
    cexhz_view = 0;

    auto khy{dz * std::atan(dy / dy) / dy};
    auto khz{dy * std::atan(dz / dy) / dz};
    auto chyex_view{xt::view(calculationParamPtr()->fdtdCoefficient()->chyex(),
                             x_range, js, xt::range(ks - 1, ks + 1))};
    chyex_view = khy * 2 * chyex_view / std::log(dz / _radius);
    auto chzex_view{xt::view(calculationParamPtr()->fdtdCoefficient()->chzex(),
                             x_range, xt::range(js - 1, js + 1), ks)};
    chzex_view = khz * 2 * chzex_view / std::log(dy / _radius);

    return;
  }

  if (_axis == Axis::XYZ::Y) {
    auto ceye_view{xt::view(calculationParamPtr()->fdtdCoefficient()->ceye(),
                            is, y_range, ks)};
    auto ceyhx_view{xt::view(calculationParamPtr()->fdtdCoefficient()->ceyhx(),
                             is, y_range, ks)};
    auto ceyhz_view{xt::view(calculationParamPtr()->fdtdCoefficient()->ceyhz(),
                             is, y_range, ks)};
    ceye_view = 0;
    ceyhx_view = 0;
    ceyhz_view = 0;

    auto khz{dx * std::atan(dz / dx) / dz};
    auto khx{dz * std::atan(dx / dz) / dx};
    auto chzey_view{xt::view(calculationParamPtr()->fdtdCoefficient()->chzey(),
                             xt::range(is - 1, is + 1), y_range, ks)};
    chzey_view = khz * 2 * chzey_view / std::log(dx / _radius);

    auto chxey_view{xt::view(calculationParamPtr()->fdtdCoefficient()->chxey(),
                             is, y_range, xt::range(ks - 1, ks + 1))};
    chxey_view = khx * 2 * chxey_view / std::log(dz / _radius);
    return;
  }

  if (_axis == Axis::XYZ::Z) {
    auto ceze_view{xt::view(calculationParamPtr()->fdtdCoefficient()->ceze(),
                            is, js, z_range)};
    auto cezhx_view{xt::view(calculationParamPtr()->fdtdCoefficient()->cezhx(),
                             is, js, z_range)};
    auto cezhy_view{xt::view(calculationParamPtr()->fdtdCoefficient()->cezhy(),
                             is, js, z_range)};
    ceze_view = 0;
    cezhx_view = 0;
    cezhy_view = 0;

    auto khx{dy * std::atan(dx / dy) / dx};
    auto khy{dx * std::atan(dy / dx) / dy};
    auto chxey_view{xt::view(calculationParamPtr()->fdtdCoefficient()->chxey(),
                             is, xt::range(js - 1, js + 1), z_range)};
    chxey_view = khx * 2 * chxey_view / std::log(dy / _radius);

    auto chyex_view{xt::view(calculationParamPtr()->fdtdCoefficient()->chyex(),
                             xt::range(is - 1, is + 1), js, z_range)};
    chyex_view = khy * 2 * chyex_view / std::log(dx / _radius);
    return;
  }
}

bool ThinWire::isEnoughThin() const {
  auto dx{gridSpacePtr()->basedDx()};
  auto dy{gridSpacePtr()->basedDy()};
  auto dz{gridSpacePtr()->basedDz()};
  if (_axis == Axis::XYZ::X) {
    return _radius < dy / 2 && _radius < dz / 2;
  }
  if (_axis == Axis::XYZ::Y) {
    return _radius < dx / 2 && _radius < dz / 2;
  }
  if (_axis == Axis::XYZ::Z) {
    return _radius < dx / 2 && _radius < dy / 2;
  }
  return false;
}

}  // namespace xfdtd
