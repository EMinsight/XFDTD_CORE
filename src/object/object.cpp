#include <xfdtd/common/type_define.h>
#include <xfdtd/grid_space/grid_space.h>
#include <xfdtd/material/ade_method/ade_method.h>
#include <xfdtd/material/dispersive_material.h>
#include <xfdtd/object/object.h>
#include <xfdtd/shape/shape.h>

#include <algorithm>
#include <cstdlib>
#include <execution>
#include <memory>
#include <utility>

#include "corrector/corrector.h"

namespace xfdtd {

Object::Object(std::string name, std::unique_ptr<Shape> shape,
               std::shared_ptr<Material> material)
    : _name{std::move(name)},
      _shape{std::move(shape)},
      _material{std::move(material)} {}

Object::~Object() = default;

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

  _grid_box = std::make_unique<GridBox>(
      _grid_space->getGridBoxWithoutCheck(_shape.get()));
  _global_grid_box =
      _grid_space->globalGridSpace()->getGridBoxWithoutCheck(shape().get());
}

void Object::correctMaterialSpace(Index index) {
  defaultCorrectMaterialSpace(index);
}

void Object::correctUpdateCoefficient() {}

void Object::handleDispersion(
    std::shared_ptr<ADEMethodStorage> ade_method_storage) {
  auto dispersion{_material->dispersion()};
  if (!dispersion) {
    return;
  }

  auto nx = _grid_space->sizeX();
  auto ny = _grid_space->sizeY();
  auto nz = _grid_space->sizeZ();

  auto linear_dispersive_material =
      dynamic_cast<LinearDispersiveMaterial*>(_material.get());

  std::for_each(std::execution::par_unseq,
                _grid_space->gridWithMaterial().begin(),
                _grid_space->gridWithMaterial().end(),
                [linear_dispersive_material, nx, ny, nz, ade_method_storage,
                 grid_space = _grid_space,
                 calculation_param = _calculation_param, this](auto&& grid) {
                  auto material_index = grid.materialIndex();
                  if (material_index != materialIndex()) {
                    return;
                  }

                  ade_method_storage->correctCoeff(
                      grid.i(), grid.j(), grid.k(), *linear_dispersive_material,
                      grid_space, calculation_param);
                });
}

void Object::correctE() {}

void Object::correctH() {}

std::unique_ptr<Corrector> Object::generateCorrector(const Task<Index>& task) {
  return nullptr;
}

std::string Object::name() const { return _name; }

const std::unique_ptr<Shape>& Object::shape() const { return _shape; }

auto Object::setMaterialIndex(Index index) -> void { _material_index = index; }

void Object::defaultCorrectMaterialSpace(Index index) {
  setMaterialIndex(index);
  auto em_property{_material->emProperty()};
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

  // // remove const
  auto g_variety = std::const_pointer_cast<GridSpace>(_grid_space);

  auto nx = _grid_space->sizeX();
  auto ny = _grid_space->sizeY();
  auto nz = _grid_space->sizeZ();

  std::for_each(std::execution::par_unseq,
                g_variety->gridWithMaterial().begin(),
                g_variety->gridWithMaterial().end(),
                [index, &grid_space = _grid_space, &shape = _shape, &eps_x,
                 &eps_y, &eps_z, &mu_x, &mu_y, &mu_z, &sigma_e_x, &sigma_e_y,
                 &sigma_e_z, &sigma_m_x, &sigma_m_y, &sigma_m_z, eps, mu,
                 sigma_e, sigma_m](auto&& g) {
                  if (!shape->isInside(grid_space->getGridCenterVector(g),
                                       grid_space->eps())) {
                    return;
                  }
                  g.setMaterialIndex(index);
                  auto i = g.i();
                  auto j = g.j();
                  auto k = g.k();
                  eps_x(i, j, k) = eps;
                  eps_y(i, j, k) = eps;
                  eps_z(i, j, k) = eps;
                  mu_x(i, j, k) = mu;
                  mu_y(i, j, k) = mu;
                  mu_z(i, j, k) = mu;
                  sigma_e_x(i, j, k) = sigma_e;
                  sigma_e_y(i, j, k) = sigma_e;
                  sigma_e_z(i, j, k) = sigma_e;
                  sigma_m_x(i, j, k) = sigma_m;
                  sigma_m_y(i, j, k) = sigma_m;
                  sigma_m_z(i, j, k) = sigma_m;
                });
}

auto Object::materialIndex() const -> Index { return _material_index; }

Shape* Object::shapePtr() { return _shape.get(); }

Material* Object::materialPtr() { return _material.get(); }

void Object::initTimeDependentVariable() {}

const GridSpace* Object::gridSpacePtr() const { return _grid_space.get(); }

CalculationParam* Object::calculationParamPtr() {
  return _calculation_param.get();
}

EMF* Object::emfPtr() { return _emf.get(); }

GridBox* Object::gridBoxPtr() const { return _grid_box.get(); }

GridBox Object::globalGridBox() const { return _global_grid_box; }

}  // namespace xfdtd
