#include <xfdtd/common/type_define.h>
#include <xfdtd/grid_space/grid_space.h>
#include <xfdtd/material/dispersive_material.h>
#include <xfdtd/object/object.h>
#include <xfdtd/shape/shape.h>

#include <algorithm>
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
  _global_grid_box = _grid_space->globalGridSpace()->getGridBox(_shape.get());
}

void Object::correctMaterialSpace(std::size_t index) {
  defaultCorrectMaterialSpace(index);
}

void Object::correctUpdateCoefficient() {
  auto dispersion{_material->dispersion()};
  if (dispersion) {
    handleDispersion();
    return;
  }
}

void Object::correctE() {}

void Object::correctH() {}

std::unique_ptr<Corrector> Object::generateCorrector(
    const Task<std::size_t>& task) {
  return nullptr;
}

std::string Object::name() const { return _name; }

const std::unique_ptr<Shape>& Object::shape() const { return _shape; }

void Object::defaultCorrectMaterialSpace(std::size_t index) {
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
          if (!shape->isInside(grid_space->getGridCenterVector({i, j, k}))) {
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
                  if (!shape->isInside(grid_space->getGridCenterVector(*g))) {
                    return;
                  }
                  g->setMaterialIndex(index);
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

Shape* Object::shapePtr() { return _shape.get(); }

Material* Object::materialPtr() { return _material.get(); }

void Object::handleDispersion() {
  if (_grid_space->type() != GridSpace::Type::UNIFORM) {
    throw XFDTDObjectException{
        "handleDispersion(): Non-uniform grid space is not supported yet"};
  }

  auto mask{_grid_space->getShapeMask(_shape.get())};
  auto& cexe = _calculation_param->fdtdCoefficient()->cexe();
  auto& cexhy = _calculation_param->fdtdCoefficient()->cexhy();
  auto& cexhz = _calculation_param->fdtdCoefficient()->cexhz();
  auto& ceye = _calculation_param->fdtdCoefficient()->ceye();
  auto& ceyhz = _calculation_param->fdtdCoefficient()->ceyhz();
  auto& ceyhx = _calculation_param->fdtdCoefficient()->ceyhx();
  auto& ceze = _calculation_param->fdtdCoefficient()->ceze();
  auto& cezhx = _calculation_param->fdtdCoefficient()->cezhx();
  auto& cezhy = _calculation_param->fdtdCoefficient()->cezhy();

  auto correct_func = [&grid_space = _grid_space, &shape = _shape](
                          Index is, Index ie, Index js, Index je, Index ks,
                          Index ke, const auto& value, auto&& arr) {
    for (auto i{is}; i < ie; ++i) {
      for (auto j{js}; j < je; ++j) {
        for (auto k{ks}; k < ke; ++k) {
          if (!shape->isInside(grid_space->getGridCenterVector({i, j, k}))) {
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

  // Support for ADE-FDTD
  if (auto lorentz_medium = dynamic_cast<LorentzMedium*>(_material.get());
      lorentz_medium != nullptr) {
    lorentz_medium->calculateCoeff(gridSpacePtr(), calculationParamPtr(),
                                   emfPtr());

    const auto& c2 = lorentz_medium->coeffForADE()._c2;
    const auto& c3 = lorentz_medium->coeffForADE()._c3;
    const auto& dx = _grid_space->basedDx();
    const auto& dy = _grid_space->basedDy();
    const auto& dz = _grid_space->basedDz();

    correct_func(0, nx, 0, ny, 0, nz, c2, cexe);
    correct_func(0, nx, 0, ny, 0, nz, -c3 / dz, cexhy);
    correct_func(0, nx, 0, ny, 0, nz, c3 / dy, cexhz);
    correct_func(0, nx, 0, ny, 0, nz, c2, ceye);
    correct_func(0, nx, 0, ny, 0, nz, -c3 / dx, ceyhz);
    correct_func(0, nx, 0, ny, 0, nz, c3 / dz, ceyhx);
    correct_func(0, nx, 0, ny, 0, nz, c2, ceze);
    correct_func(0, nx, 0, ny, 0, nz, -c3 / dy, cezhx);
    correct_func(0, nx, 0, ny, 0, nz, c3 / dx, cezhy);
    return;
  }

  if (auto drude_medium = dynamic_cast<DrudeMedium*>(_material.get());
      drude_medium != nullptr) {
    drude_medium->calculateCoeff(gridSpacePtr(), calculationParamPtr(),
                                 emfPtr());

    const auto& a = drude_medium->coeffForADE()._a;
    const auto& b = drude_medium->coeffForADE()._b;
    const auto& dx = _grid_space->basedDx();
    const auto& dy = _grid_space->basedDy();
    const auto& dz = _grid_space->basedDz();

    correct_func(0, nx, 0, ny, 0, nz, a, cexe);
    correct_func(0, nx, 0, ny, 0, nz, -b / dz, cexhy);
    correct_func(0, nx, 0, ny, 0, nz, b / dy, cexhz);
    correct_func(0, nx, 0, ny, 0, nz, a, ceye);
    correct_func(0, nx, 0, ny, 0, nz, -b / dx, ceyhz);
    correct_func(0, nx, 0, ny, 0, nz, b / dz, ceyhx);
    correct_func(0, nx, 0, ny, 0, nz, a, ceze);
    correct_func(0, nx, 0, ny, 0, nz, -b / dy, cezhx);
    correct_func(0, nx, 0, ny, 0, nz, b / dx, cezhy);
    return;
  }

  if (auto debye_medium = dynamic_cast<DebyeMedium*>(_material.get());
      debye_medium != nullptr) {
    debye_medium->calculateCoeff(gridSpacePtr(), calculationParamPtr(),
                                 emfPtr());

    const auto& a = debye_medium->coeffForADE()._a;
    const auto& b = debye_medium->coeffForADE()._b;
    const auto& dx = _grid_space->basedDx();
    const auto& dy = _grid_space->basedDy();
    const auto& dz = _grid_space->basedDz();

    correct_func(0, nx, 0, ny, 0, nz, a, cexe);
    correct_func(0, nx, 0, ny, 0, nz, -b / dz, cexhy);
    correct_func(0, nx, 0, ny, 0, nz, b / dy, cexhz);
    correct_func(0, nx, 0, ny, 0, nz, a, ceye);
    correct_func(0, nx, 0, ny, 0, nz, -b / dx, ceyhz);
    correct_func(0, nx, 0, ny, 0, nz, b / dz, ceyhx);
    correct_func(0, nx, 0, ny, 0, nz, a, ceze);
    correct_func(0, nx, 0, ny, 0, nz, -b / dy, cezhx);
    correct_func(0, nx, 0, ny, 0, nz, b / dx, cezhy);
    return;
  }
}

void Object::initTimeDependentVariable() {}

const GridSpace* Object::gridSpacePtr() const { return _grid_space.get(); }

CalculationParam* Object::calculationParamPtr() {
  return _calculation_param.get();
}

EMF* Object::emfPtr() { return _emf.get(); }

GridBox* Object::gridBoxPtr() const { return _grid_box.get(); }

GridBox Object::globalGridBox() const { return _global_grid_box; }

}  // namespace xfdtd
