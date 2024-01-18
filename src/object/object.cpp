#include <xfdtd/object/object.h>

#include <memory>
#include <utility>
#include <xtensor.hpp>

#include "xfdtd/shape/cube.h"
#include "xfdtd/shape/shape.h"

namespace xfdtd {

Object::Object(std::string name, std::unique_ptr<Shape> shape,
               std::unique_ptr<Material> material)
    : _name{std::move(name)},
      _shape{std::move(shape)},
      _material{std::move(material)} {}

std::string Object::toString() const {
  return "Object{name: " + _name + ", shape: " + _shape->toString() +
         ", material: " + _material->toString() + "}";
}

void Object::init(std::shared_ptr<const GridSpace> grid_space,
                  std::shared_ptr<CalculationParam> calculation_param,
                  std::shared_ptr<EMF> emf) {
  _grid_space = std::move(grid_space);
  _calculation_param = std::move(calculation_param);
  _emf = std::move(emf);

  _grid_box = std::make_unique<GridBox>(_grid_space->getGridBox(_shape.get()));
}

void Object::correctMaterialSpace() { defaultCorrectMaterialSpace(); }

void Object::correctUpdateCoefficient() {}

void Object::correctE() {}

void Object::correctH() {}

std::string Object::name() const { return _name; }

const std::unique_ptr<Shape>& Object::shape() const { return _shape; }

const std::unique_ptr<Material>& Object::material() const { return _material; }

void Object::defaultCorrectMaterialSpace() {
  auto dispersion{_material->dispersion()};
  if (dispersion) {
    handleDispersion();
    return;
  }

  auto em_property{_material->emProperty()};
  auto shape{_shape.get()};
  auto mask{_grid_space->getShapeMask(_shape.get())};
  auto eps{em_property.epsilon()};
  auto mu{em_property.mu()};
  auto sigma_e{em_property.sigmaE()};
  auto sigma_m{em_property.sigmaM()};
  auto& eps_x{_calculation_param->materialParam()->epsX()};
  auto& eps_y{_calculation_param->materialParam()->epsY()};
  auto& eps_z{_calculation_param->materialParam()->epsZ()};
  auto& mu_x{_calculation_param->materialParam()->muX()};
  auto& mu_y{_calculation_param->materialParam()->muY()};
  auto& mu_z{_calculation_param->materialParam()->muZ()};
  auto& sigma_e_x{_calculation_param->materialParam()->sigmaEX()};
  auto& sigma_e_y{_calculation_param->materialParam()->sigmaEY()};
  auto& sigma_e_z{_calculation_param->materialParam()->sigmaEZ()};
  auto& sigma_m_x{_calculation_param->materialParam()->sigmaMX()};
  auto& sigma_m_y{_calculation_param->materialParam()->sigmaMY()};
  auto& sigma_m_z{_calculation_param->materialParam()->sigmaMZ()};

  if (dynamic_cast<Cube*>(shape) != nullptr) {
    auto box{_grid_space->getGridBox(shape)};
    auto is{box.origin().i()};
    auto js{box.origin().j()};
    auto ks{box.origin().k()};
    auto ie{box.end().i()};
    auto je{box.end().j()};
    auto ke{box.end().k()};
    xt::view(eps_x, xt::range(is, ie), xt::range(js, je), xt::range(ks, ke)) =
        eps;
    xt::view(eps_y, xt::range(is, ie), xt::range(js, je), xt::range(ks, ke)) =
        eps;
    xt::view(eps_z, xt::range(is, ie), xt::range(js, je), xt::range(ks, ke)) =
        eps;

    xt::view(mu_x, xt::range(is, ie), xt::range(js, je), xt::range(ks, ke)) =
        mu;
    xt::view(mu_y, xt::range(is, ie), xt::range(js, je), xt::range(ks, ke)) =
        mu;
    xt::view(mu_z, xt::range(is, ie), xt::range(js, je), xt::range(ks, ke)) =
        mu;
    xt::view(sigma_e_x, xt::range(is, ie), xt::range(js, je),
             xt::range(ks, ke)) = sigma_e;
    xt::view(sigma_e_y, xt::range(is, ie), xt::range(js, je),
             xt::range(ks, ke)) = sigma_e;
    xt::view(sigma_e_z, xt::range(is, ie), xt::range(js, je),
             xt::range(ks, ke)) = sigma_e;
    xt::view(sigma_m_x, xt::range(is, ie), xt::range(js, je),
             xt::range(ks, ke)) = sigma_m;
    xt::view(sigma_m_y, xt::range(is, ie), xt::range(js, je),
             xt::range(ks, ke)) = sigma_m;
    xt::view(sigma_m_z, xt::range(is, ie), xt::range(js, je),
             xt::range(ks, ke)) = sigma_m;
    return;
  }

  xt::filter(eps_x, mask) = eps;
  xt::filter(eps_y, mask) = eps;
  xt::filter(eps_z, mask) = eps;
  xt::filter(mu_x, mask) = mu;
  xt::filter(mu_y, mask) = mu;
  xt::filter(mu_z, mask) = mu;
  xt::filter(sigma_e_x, mask) = sigma_e;
  xt::filter(sigma_e_y, mask) = sigma_e;
  xt::filter(sigma_e_z, mask) = sigma_e;
  xt::filter(sigma_m_x, mask) = sigma_m;
  xt::filter(sigma_m_y, mask) = sigma_m;
  xt::filter(sigma_m_z, mask) = sigma_m;
}

Shape* Object::shapePtr() { return _shape.get(); }

Material* Object::materialPtr() { return _material.get(); }

void Object::handleDispersion() {
  throw XFDTDObjectException{"Dispersion is not supported yet"};
}

const GridSpace* Object::gridSpacePtr() const { return _grid_space.get(); }

CalculationParam* Object::calculationParamPtr() {
  return _calculation_param.get();
}

EMF* Object::emfPtr() { return _emf.get(); }

GridBox* Object::gridBoxPtr() const { return _grid_box.get(); }

}  // namespace xfdtd
