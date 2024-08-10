#include <xfdtd/common/type_define.h>
#include <xfdtd/grid_space/grid_space.h>
#include <xfdtd/material/ade_method/ade_method.h>
#include <xfdtd/material/dispersive_material.h>
#include <xfdtd/object/object.h>
#include <xfdtd/shape/shape.h>

#include <algorithm>
#include <cstdlib>
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

  auto num_pole = linear_dispersive_material;

  for (Index i{0}; i < nx; ++i) {
    for (Index j{0}; j < ny; ++j) {
      for (Index k{0}; k < nz; ++k) {
        if (!_shape->isInside(_grid_space->getGridCenterVector({i, j, k}),
                              _grid_space->eps())) {
          continue;
        }

        const auto& grid = _grid_space->gridWithMaterial()(i, j, k);
        auto material_index = grid.materialIndex();
        if (material_index != materialIndex()) {
          continue;
        }

        ade_method_storage->correctCoeff(i, j, k, *linear_dispersive_material,
                                         _grid_space, _calculation_param);
      }
    }
  }
}

void Object::correctE() {}

void Object::correctH() {}

std::unique_ptr<Corrector> Object::generateCorrector(const Task<Index>& task) {
  return nullptr;
}

std::string Object::name() const { return _name; }

const std::unique_ptr<Shape>& Object::shape() const { return _shape; }

void Object::defaultCorrectMaterialSpace(Index index) {
  _material_index = index;
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

  auto correct_func = [&grid_space = _grid_space, &shape = _shape](
                          Index is, Index ie, Index js, Index je, Index ks,
                          Index ke, const auto& value, auto&& arr) {
    for (auto i{is}; i < ie; ++i) {
      for (auto j{js}; j < je; ++j) {
        for (auto k{ks}; k < ke; ++k) {
          if (!shape->isInside(grid_space->getGridCenterVector({i, j, k}),
                               grid_space->eps())) {
            continue;
          }
          arr(i, j, k) = value;
        }
      }
    }
  };

  auto nx = _grid_space->sizeX();
  auto ny = _grid_space->sizeY();
  auto nz = _grid_space->sizeZ();

  std::for_each(g_variety->gridWithMaterial().begin(),
                g_variety->gridWithMaterial().end(),
                [index, &grid_space = _grid_space, &shape = _shape](auto&& g) {
                  if (!shape->isInside(grid_space->getGridCenterVector(g),
                                       grid_space->eps())) {
                    return;
                  }
                  g.setMaterialIndex(index);
                });
  correct_func(0, nx, 0, ny, 0, nz, eps, eps_x);
  correct_func(0, nx, 0, ny, 0, nz, eps, eps_y);
  correct_func(0, nx, 0, ny, 0, nz, eps, eps_z);
  correct_func(0, nx, 0, ny, 0, nz, mu, mu_x);
  correct_func(0, nx, 0, ny, 0, nz, mu, mu_y);
  correct_func(0, nx, 0, ny, 0, nz, mu, mu_z);
  correct_func(0, nx, 0, ny, 0, nz, sigma_e, sigma_e_x);
  correct_func(0, nx, 0, ny, 0, nz, sigma_e, sigma_e_y);
  correct_func(0, nx, 0, ny, 0, nz, sigma_e, sigma_e_z);
  correct_func(0, nx, 0, ny, 0, nz, sigma_m, sigma_m_x);
  correct_func(0, nx, 0, ny, 0, nz, sigma_m, sigma_m_y);
  correct_func(0, nx, 0, ny, 0, nz, sigma_m, sigma_m_z);
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
