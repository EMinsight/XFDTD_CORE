#include <xfdtd/boundary/pml.h>
#include <xfdtd/calculation_param/calculation_param.h>
#include <xfdtd/coordinate_system/coordinate_system.h>
#include <xfdtd/electromagnetic_field/electromagnetic_field.h>
#include <xfdtd/grid_space/grid_space.h>
#include <xfdtd/common/constant.h>
#include <xfdtd/util/transform.h>

#include <cassert>
#include <cstddef>
#include <cstdlib>
#include <memory>
#include <string>
#include <tuple>
#include <utility>
#include <xtensor.hpp>

#include "boundary/pml_corrector.h"
#include "corrector/corrector.h"

namespace xfdtd {

PML::PML(int thickness, Axis::Direction direction, int order,
         Real sigma_ratio, Real alpha_min, Real alpha_max,
         Real kappa_max)
    : _thickness{thickness},
      _n{static_cast<std::size_t>(std::abs(thickness))},
      _direction{direction},
      _main_axis{Axis::fromDirectionToXYZ(direction)},
      _order{order},
      _sigma_ratio{sigma_ratio},
      _alpha_min{alpha_min},
      _alpha_max{alpha_max},
      _kappa_max{kappa_max} {}

int PML::thickness() const { return _thickness; }

Axis::Direction PML::direction() const { return _direction; }

Axis::XYZ PML::subAxisA() const {
  switch (_direction) {
    case Axis::Direction::XP:
    case Axis::Direction::XN:
      return Axis::XYZ::Y;
    case Axis::Direction::YP:
    case Axis::Direction::YN:
      return Axis::XYZ::Z;
    case Axis::Direction::ZP:
    case Axis::Direction::ZN:
      return Axis::XYZ::X;
    default:
      throw XFDTDPMLException("Invalid direction");
  }
}

Axis::XYZ PML::subAxisB() const {
  switch (_direction) {
    case Axis::Direction::XP:
    case Axis::Direction::XN:
      return Axis::XYZ::Z;
    case Axis::Direction::YP:
    case Axis::Direction::YN:
      return Axis::XYZ::X;
    case Axis::Direction::ZP:
    case Axis::Direction::ZN:
      return Axis::XYZ::Y;
    default:
      throw XFDTDPMLException("Invalid direction");
  }
}

Axis::XYZ PML::mainAxis() const { return _main_axis; };

std::size_t PML::globalENodeStartIndexMainAxis() const {
  return _global_e_start_index;
}

std::size_t PML::globalHNodeStartIndexMainAxis() const {
  return _global_h_start_index;
}

std::size_t PML::nodeENodeStartIndexMainAxis() const {
  return _node_e_start_index;
}

std::size_t PML::nodeHNodeStartIndexMainAxis() const {
  return _node_h_start_index;
}

std::size_t PML::n() const { return _n; }

std::size_t PML::nodeN() const { return _pml_node_task_abc.zRange().size(); }

const Array1D<Real>& PML::globalESize() const { return _global_e_size; }

const Array1D<Real>& PML::globalHSize() const { return _global_h_size; }

Array3D<Real>& PML::eaF() {
  switch (mainAxis()) {
    case Axis::XYZ::X:
      return emfPtr()->ey();
    case Axis::XYZ::Y:
      return emfPtr()->ez();
    case Axis::XYZ::Z:
      return emfPtr()->ex();
    default:
      throw XFDTDPMLException("Invalid direction");
  }
}

Array3D<Real>& PML::ebF() {
  switch (mainAxis()) {
    case Axis::XYZ::X:
      return emfPtr()->ez();
    case Axis::XYZ::Y:
      return emfPtr()->ex();
    case Axis::XYZ::Z:
      return emfPtr()->ey();
    default:
      throw XFDTDPMLException("Invalid direction");
  }
}

Array3D<Real>& PML::haF() {
  switch (Axis::fromDirectionToXYZ(_direction)) {
    case Axis::XYZ::X:
      return emfPtr()->hy();
    case Axis::XYZ::Y:
      return emfPtr()->hz();
    case Axis::XYZ::Z:
      return emfPtr()->hx();
    default:
      throw XFDTDPMLException("Invalid direction");
  }
}

Array3D<Real>& PML::hbF() {
  switch (mainAxis()) {
    case Axis::XYZ::X:
      return emfPtr()->hz();
    case Axis::XYZ::Y:
      return emfPtr()->hx();
    case Axis::XYZ::Z:
      return emfPtr()->hy();
    default:
      throw XFDTDPMLException("Invalid direction");
  }
}

void PML::init(std::shared_ptr<const GridSpace> grid_space,
               std::shared_ptr<CalculationParam> calculation_param,
               std::shared_ptr<EMF> emf) {
  defaultInit(std::move(grid_space), std::move(calculation_param),
              std::move(emf));

  auto global_box{gridSpacePtr()->globalGridSpace()->box()};
  switch (_direction) {
    case Axis::Direction::XN: {
      _global_e_start_index = global_box.origin().i() + 1;
      _global_h_start_index = global_box.origin().i();
      break;
    }
    case Axis::Direction::XP: {
      _global_e_start_index = global_box.end().i() - n();
      _global_h_start_index = global_box.end().i() - n();
      break;
    }
    case Axis::Direction::YN: {
      _global_e_start_index = global_box.origin().j() + 1;
      _global_h_start_index = global_box.origin().j();
      break;
    }
    case Axis::Direction::YP: {
      _global_e_start_index = global_box.end().j() - n();
      _global_h_start_index = global_box.end().j() - n();
      break;
    }
    case Axis::Direction::ZN: {
      _global_e_start_index = global_box.origin().k() + 1;
      _global_h_start_index = global_box.origin().k();
      break;
    }
    case Axis::Direction::ZP: {
      _global_e_start_index = global_box.end().k() - n();
      _global_h_start_index = global_box.end().k() - n();
      break;
    }
    default:
      throw XFDTDPMLException("Invalid direction");
  }

  // _na = gridSpacePtr()->eSize(subAxisA()).size();
  // _nb = gridSpacePtr()->eSize(subAxisB()).size();
  const std::shared_ptr<const GridSpace> global_space =
      gridSpacePtr()->globalGridSpace();
  _global_na = global_space->eSize(subAxisA()).size();
  _global_nb = global_space->eSize(subAxisB()).size();

  auto node_global_box = gridSpacePtr()->globalBox();
  auto node_global_task = makeIndexTask(
      makeIndexRange(node_global_box.origin().i(),
                              node_global_box.end().i()),
      makeIndexRange(node_global_box.origin().j(),
                              node_global_box.end().j()),
      makeIndexRange(node_global_box.origin().k(),
                              node_global_box.end().k()));

  _pml_global_task_abc = makeIndexTask(
      makeIndexRange(0, _global_na),
      makeIndexRange(0, _global_nb),
      makeIndexRange(globalHNodeStartIndexMainAxis(),
                              globalHNodeStartIndexMainAxis() + n()));
  auto [x, y, z] =
      transform::aBCToXYZ(std::make_tuple(_pml_global_task_abc.xRange(),
                                          _pml_global_task_abc.yRange(),
                                          _pml_global_task_abc.zRange()),
                          mainAxis());
  _pml_global_task = makeIndexTask(x, y, z);

  auto t = taskIntersection(node_global_task, _pml_global_task);
  if (!t.has_value()) {
    _pml_node_task = makeIndexTask(makeIndexRange(0, 0),
                                            makeIndexRange(0, 0),
                                            makeIndexRange(0, 0));
    _pml_node_task_abc = makeIndexTask(makeIndexRange(0, 0),
                                                makeIndexRange(0, 0),
                                                makeIndexRange(0, 0));
    return;
  }

  _pml_node_task = t.value();
  _pml_node_task = makeIndexTask(
      _pml_node_task.xRange() - node_global_box.origin().i(),
      _pml_node_task.yRange() - node_global_box.origin().j(),
      _pml_global_task.zRange() - node_global_box.origin().k());

  auto [a, b, c] = transform::xYZToABC(
      std::make_tuple(_pml_node_task.xRange(), _pml_node_task.yRange(),
                      _pml_node_task.zRange()),
      mainAxis());
  _pml_node_task_abc = makeIndexTask(a, b, c);

  if (Axis::directionNegative(direction())) {
    _node_e_start_index = _pml_node_task_abc.zRange().start() + 1;
    _node_h_start_index = _pml_node_task_abc.zRange().start();
  } else {
    _node_e_start_index = _pml_node_task_abc.zRange().start();
    _node_h_start_index = _pml_node_task_abc.zRange().start();
  }
}

void PML::correctMaterialSpace() {}

void PML::correctUpdateCoefficient() {
  if (!_pml_node_task.valid()) {
    return;
  }

  std::shared_ptr<const GridSpace> global_space =
      gridSpacePtr()->globalGridSpace();

  _global_h_size = xt::view(global_space->hSize(mainAxis()),
                            xt::range(globalENodeStartIndexMainAxis(),
                                      globalENodeStartIndexMainAxis() + n()));
  _global_e_size = xt::view(global_space->eSize(mainAxis()),
                            xt::range(globalHNodeStartIndexMainAxis(),
                                      globalHNodeStartIndexMainAxis() + n()));

  calRecursiveConvolutionCoeff();

  switch (mainAxis()) {
    case Axis::XYZ::X:
      correctCoefficientX();
      return;
    case Axis::XYZ::Y:
      correctCoefficientY();
      return;
    case Axis::XYZ::Z:
      correctCoefficientZ();
      return;
    default:
      throw XFDTDPMLException("Invalid direction");
  }
}

std::unique_ptr<Corrector> PML::generateDomainCorrector(
    const Task<std::size_t>& task) {
  if (!_pml_node_task.valid() || !intersected(task, _pml_node_task)) {
    return nullptr;
  }

  auto t = taskIntersection(task, _pml_node_task);
  if (!t.has_value()) {
    return nullptr;
  }

  auto [a, b, c_h_n] = transform::xYZToABC(
      std::make_tuple(t.value().xRange(), t.value().yRange(),
                      t.value().zRange()),
      mainAxis());
  auto t_abc = makeIndexTask(a, b, c_h_n);

  auto a_s = t_abc.xRange().start();
  auto a_n = t_abc.xRange().size();
  auto b_s = t_abc.yRange().start();
  auto b_n = t_abc.yRange().size();
  auto c_s = t_abc.zRange().start();
  auto c_n = t_abc.zRange().size();
  auto c_e_s = (Axis::directionNegative(direction())) ? (c_s + 1) : (c_s);

  auto offset = gridSpacePtr()->globalBox().origin();

  switch (mainAxis()) {
    case Axis::XYZ::X:
      return std::make_unique<PMLCorrectorX>(
          globalENodeStartIndexMainAxis(), globalHNodeStartIndexMainAxis(),
          nodeENodeStartIndexMainAxis(), nodeHNodeStartIndexMainAxis(), a_s,
          a_n, b_s, b_n, c_e_s, c_n, c_s, c_n, offset.i(), _coeff_a_e,
          _coeff_b_e, _coeff_a_h, _coeff_b_h, _c_ea_psi_hb, _c_eb_psi_ha,
          _c_ha_psi_eb, _c_hb_psi_ea, _ea_psi_hb, _eb_psi_ha, _ha_psi_eb,
          _hb_psi_ea, eaF(), ebF(), haF(), hbF());
    case Axis::XYZ::Y:
      return std::make_unique<PMLCorrectorY>(
          globalENodeStartIndexMainAxis(), globalHNodeStartIndexMainAxis(),
          nodeENodeStartIndexMainAxis(), nodeHNodeStartIndexMainAxis(), a_s,
          a_n, b_s, b_n, c_e_s, c_n, c_s, c_n, offset.j(), _coeff_a_e,
          _coeff_b_e, _coeff_a_h, _coeff_b_h, _c_ea_psi_hb, _c_eb_psi_ha,
          _c_ha_psi_eb, _c_hb_psi_ea, _ea_psi_hb, _eb_psi_ha, _ha_psi_eb,
          _hb_psi_ea, eaF(), ebF(), haF(), hbF());
    case Axis::XYZ::Z:
      return std::make_unique<PMLCorrectorZ>(
          globalENodeStartIndexMainAxis(), globalHNodeStartIndexMainAxis(),
          nodeENodeStartIndexMainAxis(), nodeHNodeStartIndexMainAxis(), a_s,
          a_n, b_s, b_n, c_e_s, c_n, c_s, c_n, offset.k(), _coeff_a_e,
          _coeff_b_e, _coeff_a_h, _coeff_b_h, _c_ea_psi_hb, _c_eb_psi_ha,
          _c_ha_psi_eb, _c_hb_psi_ea, _ea_psi_hb, _eb_psi_ha, _ha_psi_eb,
          _hb_psi_ea, eaF(), ebF(), haF(), hbF());
    default:
      throw XFDTDPMLException("Invalid direction");
  }
}
// }

void PML::correctCoefficientX() {
  const auto node_nx = _pml_node_task.xRange().size();
  const auto node_ny = _pml_node_task.yRange().size();
  const auto node_nz = _pml_node_task.zRange().size();

  const auto global_e_start = globalENodeStartIndexMainAxis();
  const auto global_h_start = globalHNodeStartIndexMainAxis();
  const auto node_e_start = nodeENodeStartIndexMainAxis();
  const auto node_h_start = nodeHNodeStartIndexMainAxis();
  const auto offset = gridSpacePtr()->globalBox().origin();

  _c_ea_psi_hb = xt::zeros<Real>({node_nx, node_ny, node_nz + 1});
  _ea_psi_hb = xt::zeros<Real>({node_nx, node_ny, node_nz + 1});

  auto& ceahb{calculationParamPtr()->fdtdCoefficient()->ceyhz()};

  _c_eb_psi_ha = xt::zeros<Real>({node_nx, node_ny + 1, node_nz});
  _eb_psi_ha = xt::zeros<Real>({node_nx, node_ny + 1, node_nz});

  auto& cebha{calculationParamPtr()->fdtdCoefficient()->cezhy()};

  for (std::size_t i{0}; i < node_nx; ++i) {
    auto node_i = i + node_e_start;
    auto global_i = node_i + offset.i();

    for (std::size_t j{0}; j < node_ny; ++j) {
      for (std::size_t k{0}; k < node_nz + 1; ++k) {
        _c_ea_psi_hb(i, j, k) =
            ceahb(node_i, j, k) * _global_h_size(global_i - global_e_start);
        ceahb(node_i, j, k) =
            ceahb(node_i, j, k) / _kappa_e(global_i - global_e_start);
      }
    }

    for (std::size_t j{0}; j < node_ny + 1; ++j) {
      for (std::size_t k{0}; k < node_nz; ++k) {
        _c_eb_psi_ha(i, j, k) =
            cebha(node_i, j, k) * _global_h_size(global_i - global_e_start);
        cebha(node_i, j, k) =
            cebha(node_i, j, k) / _kappa_e(global_i - global_e_start);
      }
    }
  }

  _c_ha_psi_eb = xt::zeros<Real>({node_nx, node_ny + 1, node_nz});
  _ha_psi_eb = xt::zeros<Real>({node_nx, node_ny + 1, node_nz});

  auto& chaeb{calculationParamPtr()->fdtdCoefficient()->chyez()};

  _c_hb_psi_ea = xt::zeros<Real>({node_nx, node_ny, node_nz + 1});
  _hb_psi_ea = xt::zeros<Real>({node_nx, node_ny, node_nz + 1});

  auto& chbea{calculationParamPtr()->fdtdCoefficient()->chzey()};

  for (std::size_t i{0}; i < node_nx; ++i) {
    auto node_i = i + node_h_start;
    auto global_i = node_i + offset.i();

    for (std::size_t j{0}; j < node_ny + 1; ++j) {
      for (std::size_t k{0}; k < node_nz; ++k) {
        _c_ha_psi_eb(i, j, k) =
            chaeb(node_i, j, k) * _global_e_size(global_i - global_h_start);
        chaeb(node_i, j, k) =
            chaeb(node_i, j, k) / _kappa_h(global_i - global_h_start);
      }
    }

    for (std::size_t j{0}; j < node_ny; ++j) {
      for (std::size_t k{0}; k < node_nz + 1; ++k) {
        _c_hb_psi_ea(i, j, k) =
            chbea(node_i, j, k) * _global_e_size(global_i - global_h_start);
        chbea(node_i, j, k) =
            chbea(node_i, j, k) / _kappa_h(global_i - global_h_start);
      }
    }
  }
}

void PML::correctCoefficientY() {
  const auto node_nx = _pml_node_task.xRange().size();
  const auto node_ny = _pml_node_task.yRange().size();
  const auto node_nz = _pml_node_task.zRange().size();

  const auto global_e_start = globalENodeStartIndexMainAxis();
  const auto global_h_start = globalHNodeStartIndexMainAxis();
  const auto node_e_start = nodeENodeStartIndexMainAxis();
  const auto node_h_start = nodeHNodeStartIndexMainAxis();
  const auto offset = gridSpacePtr()->globalBox().origin();

  _c_ea_psi_hb = xt::zeros<Real>({node_nx + 1, node_ny, node_nz});
  _ea_psi_hb = xt::zeros<Real>({node_nx + 1, node_ny, node_nz});
  auto& ceahb{calculationParamPtr()->fdtdCoefficient()->cezhx()};

  _c_eb_psi_ha = xt::zeros<Real>({node_nx, node_ny, node_nz + 1});
  _eb_psi_ha = xt::zeros<Real>({node_nx, node_ny, node_nz + 1});
  auto& cebha{calculationParamPtr()->fdtdCoefficient()->cexhz()};

  for (std::size_t i{0}; i < node_nx + 1; ++i) {
    for (std::size_t j{0}; j < node_ny; ++j) {
      auto node_j = j + node_e_start;
      auto global_j = node_j + offset.j();
      for (std::size_t k{0}; k < node_nz; ++k) {
        _c_ea_psi_hb(i, j, k) =
            ceahb(i, node_j, k) * _global_h_size(global_j - global_e_start);
        ceahb(i, node_j, k) =
            ceahb(i, node_j, k) / _kappa_e(global_j - global_e_start);
      }
    }
  }

  for (std::size_t i{0}; i < node_nx; ++i) {
    for (std::size_t j{0}; j < node_ny; ++j) {
      auto node_j = j + node_e_start;
      auto global_j = node_j + offset.j();
      for (std::size_t k{0}; k < node_nz + 1; ++k) {
        _c_eb_psi_ha(i, j, k) =
            cebha(i, node_j, k) * _global_h_size(global_j - global_e_start);
        cebha(i, node_j, k) =
            cebha(i, node_j, k) / _kappa_e(global_j - global_e_start);
      }
    }
  }

  _c_ha_psi_eb = xt::zeros<Real>({node_nx, node_ny, node_nz + 1});
  _ha_psi_eb = xt::zeros<Real>({node_nx, node_ny, node_nz + 1});
  auto& chaeb{calculationParamPtr()->fdtdCoefficient()->chzex()};

  _c_hb_psi_ea = xt::zeros<Real>({node_nx + 1, node_ny, node_nz});
  _hb_psi_ea = xt::zeros<Real>({node_nx + 1, node_ny, node_nz});
  auto& chbea{calculationParamPtr()->fdtdCoefficient()->chxez()};

  for (std::size_t i{0}; i < node_nx; ++i) {
    for (std::size_t j{0}; j < node_ny; ++j) {
      auto node_j = j + node_h_start;
      auto global_j = node_j + offset.j();
      for (std::size_t k{0}; k < node_nz + 1; ++k) {
        _c_ha_psi_eb(i, j, k) =
            chaeb(i, node_j, k) * _global_e_size(global_j - global_h_start);
        chaeb(i, node_j, k) =
            chaeb(i, node_j, k) / _kappa_h(global_j - global_h_start);
      }
    }
  }

  for (std::size_t i{0}; i < node_nx + 1; ++i) {
    for (std::size_t j{0}; j < node_ny; ++j) {
      auto node_j = j + node_h_start;
      auto global_j = node_j + offset.j();
      for (std::size_t k{0}; k < node_nz; ++k) {
        _c_hb_psi_ea(i, j, k) =
            chbea(i, node_j, k) * _global_e_size(global_j - global_h_start);
        chbea(i, node_j, k) =
            chbea(i, node_j, k) / _kappa_h(global_j - global_h_start);
      }
    }
  }
}

void PML::correctCoefficientZ() {
  const auto node_nx = _pml_node_task.xRange().size();
  const auto node_ny = _pml_node_task.yRange().size();
  const auto node_nz = _pml_node_task.zRange().size();

  const auto global_e_start = globalENodeStartIndexMainAxis();
  const auto global_h_start = globalHNodeStartIndexMainAxis();
  const auto node_e_start = nodeENodeStartIndexMainAxis();
  const auto node_h_start = nodeHNodeStartIndexMainAxis();
  const auto offset = gridSpacePtr()->globalBox().origin();

  _c_ea_psi_hb = xt::zeros<Real>({node_nx, node_ny + 1, node_nz});
  _ea_psi_hb = xt::zeros<Real>({node_nx, node_ny + 1, node_nz});
  auto& ceahb{calculationParamPtr()->fdtdCoefficient()->cexhy()};

  _c_eb_psi_ha = xt::zeros<Real>({node_nx + 1, node_ny, node_nz});
  _eb_psi_ha = xt::zeros<Real>({node_nx + 1, node_ny, node_nz});
  auto& cebha{calculationParamPtr()->fdtdCoefficient()->ceyhx()};

  for (std::size_t i{0}; i < node_nx; ++i) {
    for (std::size_t j{0}; j < node_ny + 1; ++j) {
      for (std::size_t k{0}; k < node_nz; ++k) {
        auto node_k = k + node_e_start;
        auto global_k = node_k + offset.k();
        _c_ea_psi_hb(i, j, k) =
            ceahb(i, j, node_k) * _global_h_size(global_k - global_e_start);
        ceahb(i, j, node_k) =
            ceahb(i, j, node_k) / _kappa_e(global_k - global_e_start);
      }
    }
  }

  for (std::size_t i{0}; i < node_nx + 1; ++i) {
    for (std::size_t j{0}; j < node_ny; ++j) {
      for (std::size_t k{0}; k < node_nz; ++k) {
        auto node_k = k + node_e_start;
        auto global_k = node_k + offset.k();
        _c_eb_psi_ha(i, j, k) =
            cebha(i, j, node_k) * _global_h_size(global_k - global_e_start);
        cebha(i, j, node_k) =
            cebha(i, j, node_k) / _kappa_e(global_k - global_e_start);
      }
    }
  }

  _c_ha_psi_eb = xt::zeros<Real>({node_nx + 1, node_ny, node_nz});
  _ha_psi_eb = xt::zeros<Real>({node_nx + 1, node_ny, node_nz});

  auto& chaeb{calculationParamPtr()->fdtdCoefficient()->chxey()};

  _c_hb_psi_ea = xt::zeros<Real>({node_nx, node_ny + 1, node_nz});
  _hb_psi_ea = xt::zeros<Real>({node_nx, node_ny + 1, node_nz});

  auto& chbea{calculationParamPtr()->fdtdCoefficient()->chyex()};

  for (std::size_t i{0}; i < node_nx + 1; ++i) {
    for (std::size_t j{0}; j < node_ny; ++j) {
      for (std::size_t k{0}; k < node_nz; ++k) {
        auto node_k = k + node_h_start;
        auto global_k = node_k + offset.k();
        _c_ha_psi_eb(i, j, k) =
            chaeb(i, j, node_k) * _global_e_size(global_k - global_h_start);
        chaeb(i, j, node_k) =
            chaeb(i, j, node_k) / _kappa_h(global_k - global_h_start);
      }
    }
  }

  for (std::size_t i{0}; i < node_nx; ++i) {
    for (std::size_t j{0}; j < node_ny + 1; ++j) {
      for (std::size_t k{0}; k < node_nz; ++k) {
        auto node_k = k + node_h_start;
        auto global_k = node_k + offset.k();
        _c_hb_psi_ea(i, j, k) =
            chbea(i, j, node_k) * _global_e_size(global_k - global_h_start);
        chbea(i, j, node_k) =
            chbea(i, j, node_k) / _kappa_h(global_k - global_h_start);
      }
    }
  }
}

void PML::calRecursiveConvolutionCoeff() {
  Real min_h_size{xt::amin(globalHSize())()};
  Real sigma_max_e = calculateSigmaMax(min_h_size);
  auto rho_e{calculateRhoE(n(), globalHSize())};
  auto sigma_e{calculateSigma(sigma_max_e, rho_e, _order)};
  auto alpha_e{calculateAlpha(_alpha_min, _alpha_max, rho_e)};

  _kappa_e = calculateKappa(_kappa_max, rho_e, _order);

  _coeff_b_e = calculateCoefficientB(sigma_e, _kappa_e, alpha_e,
                                     calculationParamPtr()->timeParam()->dt(),
                                     constant::EPSILON_0);

  _coeff_a_e = calculateCoefficientA(_coeff_b_e, sigma_e, _kappa_e, alpha_e,
                                     globalHSize());

  Real e2m{constant::MU_0 / constant::EPSILON_0};
  auto min_e_size{xt::amin(globalESize())()};
  Real sigma_max_m{calculateSigmaMax(min_e_size) * e2m};
  auto rho_m{calculateRhoM(n(), globalESize())};
  auto sigma_m{calculateSigma(sigma_max_m, rho_m, _order)};
  auto alpha_m{calculateAlpha(_alpha_min, _alpha_max, rho_m) * e2m};

  _kappa_h = {calculateKappa(_kappa_max, rho_m, _order)};

  _coeff_b_h = calculateCoefficientB(sigma_m, _kappa_h, alpha_m,
                                     calculationParamPtr()->timeParam()->dt(),
                                     constant::MU_0);

  _coeff_a_h = calculateCoefficientA(_coeff_b_h, sigma_m, _kappa_h, alpha_m,
                                     globalESize());
}

Real PML::calculateSigmaMax(Real dl) const {
  return _sigma_ratio * (_order + 1) / (150 * constant::PI * dl);
}

Array1D<Real> PML::calculateRhoE(std::size_t n,
                                      const Array1D<Real>& size) const {
  assert(n == size.size());
  auto interval{size - 0.75 * size};
  Array1D<Real> d;
  d.resize({n});
  d(0) = interval(0);
  Real sum{0};
  for (std::size_t i{1}; i < n; ++i) {
    sum += size(i - 1);
    d(i) = sum + interval(i);
  }
  sum += size(n - 1);
  if (_direction == Axis::Direction::XN || _direction == Axis::Direction::YN ||
      _direction == Axis::Direction::ZN) {
    d = xt::flip(d);
  }

  return d / sum;
}

Array1D<Real> PML::calculateRhoM(std::size_t n,
                                      const Array1D<Real>& size) const {
  assert(n == size.size());
  auto interval{size - 0.25 * size};
  Array1D<Real> d;
  d.resize({n});
  d(0) = interval(0);
  Real sum{0};
  for (std::size_t i{1}; i < n; ++i) {
    sum += size(i - 1);
    d(i) = sum + interval(i);
  }
  sum += size(n - 1);
  if (_direction == Axis::Direction::XN || _direction == Axis::Direction::YN ||
      _direction == Axis::Direction::ZN) {
    d = xt::flip(d);
  }

  return d / sum;
}

Array1D<Real> PML::calculateSigma(Real sigma_max,
                                       const Array1D<Real>& rho,
                                       std::size_t order) const {
  return sigma_max * xt::pow(rho, _order);
}

Array1D<Real> PML::calculateKappa(Real kappa_max,
                                       const Array1D<Real>& rho,
                                       std::size_t order) const {
  return 1 + (kappa_max - 1) * xt::pow(rho, _order);
}

Array1D<Real> PML::calculateAlpha(Real alpha_min, Real alpha_max,
                                       const Array1D<Real>& rho) const {
  return alpha_min + (alpha_max - alpha_min) * (1 - rho);
}

Array1D<Real> PML::calculateCoefficientA(
    const Array1D<Real>& b, const Array1D<Real>& sigma,
    const Array1D<Real>& kappa, const Array1D<Real>& alpha,
    const Array1D<Real>& dl) const {
  return (b - 1) * sigma / ((sigma + kappa * alpha) * kappa * dl);
}

Array1D<Real> PML::calculateCoefficientB(const Array1D<Real>& sigma,
                                              const Array1D<Real>& kappa,
                                              const Array1D<Real>& alpha,
                                              Real dt,
                                              Real constant) const {
  return xt::exp((-dt / constant) * (sigma / kappa + alpha));
}

Array1D<Real> PML::calculateCoeffPsi(const Array1D<Real>& coeff,
                                          const Array1D<Real>& kappa,
                                          const Array1D<Real>& dl) const {
  return coeff * dl;
}

}  // namespace xfdtd
