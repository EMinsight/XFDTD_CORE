#include <xfdtd/boundary/pml.h>

#include <cassert>
#include <cstddef>
#include <cstdlib>
#include <utility>
#include <xtensor.hpp>

#include "xfdtd/coordinate_system/coordinate_system.h"
#include "xfdtd/electromagnetic_field/electromagnetic_field.h"
#include "xfdtd/util/constant.h"

namespace xfdtd {

PML::PML(int thickness, Axis::Direction direction, int order,
         double sigma_ratio, double alpha_min, double alpha_max,
         double kappa_max)
    : _thickness{thickness},
      _n{static_cast<std::size_t>(std::abs(thickness))},
      _direction{direction},
      _main_axis{Axis::formDirectionToXYZ(direction)},
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

std::size_t PML::eNodeStartIndexMainAxis() const { return _e_start_index; }

std::size_t PML::hNodeStartIndexMainAxis() const { return _h_start_index; }

std::size_t PML::n() const { return _n; }

const xt::xarray<double>& PML::eSize() const { return _e_size; }

const xt::xarray<double>& PML::hSize() const { return _h_size; }

xt::xarray<double>& PML::eaF() {
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

xt::xarray<double>& PML::ebF() {
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

xt::xarray<double>& PML::haF() {
  switch (Axis::formDirectionToXYZ(_direction)) {
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

xt::xarray<double>& PML::hbF() {
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
  auto box{gridSpacePtr()->box()};
  switch (_direction) {
    case Axis::Direction::XN: {
      _e_start_index = box.origin().i() + 1;
      _h_start_index = box.origin().i();
      break;
    }
    case Axis::Direction::XP: {
      _e_start_index = box.end().i() - n();
      _h_start_index = box.end().i() - n();
      break;
    }
    case Axis::Direction::YN: {
      _e_start_index = box.origin().j() + 1;
      _h_start_index = box.origin().j();
      break;
    }
    case Axis::Direction::YP: {
      _e_start_index = box.end().j() - n();
      _h_start_index = box.end().j() - n();
      break;
    }
    case Axis::Direction::ZN: {
      _e_start_index = box.origin().k() + 1;
      _h_start_index = box.origin().k();
      break;
    }
    case Axis::Direction::ZP: {
      _e_start_index = box.end().k() - n();
      _h_start_index = box.end().k() - n();
      break;
    }
    default:
      throw XFDTDPMLException("Invalid direction");
  }

  _na = gridSpacePtr()->eSize(subAxisA()).size();
  _nb = gridSpacePtr()->eSize(subAxisB()).size();
}

void PML::correctMaterialSpace() {}

void PML::correctUpdateCoefficient() {
  _h_size = xt::view(
      gridSpacePtr()->hSize(mainAxis()),
      xt::range(eNodeStartIndexMainAxis(), eNodeStartIndexMainAxis() + n()));
  _e_size = xt::view(
      gridSpacePtr()->eSize(mainAxis()),
      xt::range(hNodeStartIndexMainAxis(), hNodeStartIndexMainAxis() + n()));

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

void PML::correctE() {
  const auto& hb{hbF()};
  auto& ea{eaF()};

  const auto& ha{haF()};
  auto& eb{ebF()};

  switch (mainAxis()) {
    case Axis::XYZ::X: {
      for (std::size_t i{eNodeStartIndexMainAxis()};
           i < eNodeStartIndexMainAxis() + n(); ++i) {
        auto ii{i - eNodeStartIndexMainAxis()};
        for (std::size_t j{0}; j < _na; ++j) {
          for (std::size_t k{0}; k < _nb + 1; ++k) {
            _ea_psi_hb(ii, j, k) =
                _coeff_b_e(ii) * _ea_psi_hb(ii, j, k) +
                _coeff_a_e(ii) * (hb(i, j, k) - hb(i - 1, j, k));
            ea(i, j, k) += _c_ea_psi_hb(ii, j, k) * _ea_psi_hb(ii, j, k);
          }
        }

        for (std::size_t j{0}; j < _na + 1; ++j) {
          for (std::size_t k{0}; k < _nb; ++k) {
            _eb_psi_ha(ii, j, k) =
                _coeff_b_e(ii) * _eb_psi_ha(ii, j, k) +
                _coeff_a_e(ii) * (ha(i, j, k) - ha(i - 1, j, k));
            eb(i, j, k) += _c_eb_psi_ha(ii, j, k) * _eb_psi_ha(ii, j, k);
          }
        }
      }
      return;
    }
    case Axis::XYZ::Y: {
      for (std::size_t i{0}; i < _nb + 1; ++i) {
        for (std::size_t j{eNodeStartIndexMainAxis()};
             j < eNodeStartIndexMainAxis() + n(); ++j) {
          for (std::size_t k{0}; k < _na; ++k) {
            auto jj{j - eNodeStartIndexMainAxis()};
            _ea_psi_hb(i, jj, k) =
                _coeff_b_e(jj) * _ea_psi_hb(i, jj, k) +
                _coeff_a_e(jj) * (hb(i, j, k) - hb(i, j - 1, k));
            ea(i, j, k) += _c_ea_psi_hb(i, jj, k) * _ea_psi_hb(i, jj, k);
          }
        }
      }

      for (std::size_t i{0}; i < _nb; ++i) {
        for (std::size_t j{eNodeStartIndexMainAxis()};
             j < eNodeStartIndexMainAxis() + n(); ++j) {
          for (std::size_t k{0}; k < _na + 1; ++k) {
            auto jj{j - eNodeStartIndexMainAxis()};
            _eb_psi_ha(i, jj, k) =
                _coeff_b_e(jj) * _eb_psi_ha(i, jj, k) +
                _coeff_a_e(jj) * (ha(i, j, k) - ha(i, j - 1, k));
            eb(i, j, k) += _c_eb_psi_ha(i, jj, k) * _eb_psi_ha(i, jj, k);
          }
        }
      }
      return;
    }
    case Axis::XYZ::Z: {
      for (std::size_t i{0}; i < _na; ++i) {
        for (std::size_t j{0}; j < _nb + 1; ++j) {
          for (std::size_t k{eNodeStartIndexMainAxis()};
               k < eNodeStartIndexMainAxis() + n(); ++k) {
            auto kk{k - eNodeStartIndexMainAxis()};
            _ea_psi_hb(i, j, kk) =
                _coeff_b_e(kk) * _ea_psi_hb(i, j, kk) +
                _coeff_a_e(kk) * (hb(i, j, k) - hb(i, j, k - 1));
            ea(i, j, k) += _c_ea_psi_hb(i, j, kk) * _ea_psi_hb(i, j, kk);
          }
        }
      }
      for (std::size_t i{0}; i < _na + 1; ++i) {
        for (std::size_t j{0}; j < _nb; ++j) {
          for (std::size_t k{eNodeStartIndexMainAxis()};
               k < eNodeStartIndexMainAxis() + n(); ++k) {
            auto kk{k - eNodeStartIndexMainAxis()};
            _eb_psi_ha(i, j, kk) =
                _coeff_b_e(kk) * _eb_psi_ha(i, j, kk) +
                _coeff_a_e(kk) * (ha(i, j, k) - ha(i, j, k - 1));
            eb(i, j, k) += _c_eb_psi_ha(i, j, kk) * _eb_psi_ha(i, j, kk);
          }
        }
      }
      return;
    }
    default:
      throw XFDTDPMLException("Invalid direction");
  }
}

void PML::correctH() {
  const auto& ea{eaF()};
  auto& hb{hbF()};

  const auto& eb{ebF()};
  auto& ha{haF()};

  switch (mainAxis()) {
    case Axis::XYZ::X: {
      for (std::size_t i{hNodeStartIndexMainAxis()};
           i < hNodeStartIndexMainAxis() + n(); ++i) {
        auto ii{i - hNodeStartIndexMainAxis()};
        for (std::size_t j{0}; j < _na + 1; ++j) {
          for (std::size_t k{0}; k < _nb; ++k) {
            _ha_psi_eb(ii, j, k) =
                _coeff_b_h(ii) * _ha_psi_eb(ii, j, k) +
                _coeff_a_h(ii) * (eb(i + 1, j, k) - eb(i, j, k));
            ha(i, j, k) += _c_ha_psi_eb(ii, j, k) * _ha_psi_eb(ii, j, k);
          }
        }

        for (std::size_t j{0}; j < _na; ++j) {
          for (std::size_t k{0}; k < _nb + 1; ++k) {
            _hb_psi_ea(ii, j, k) =
                _coeff_b_h(ii) * _hb_psi_ea(ii, j, k) +
                _coeff_a_h(ii) * (ea(i + 1, j, k) - ea(i, j, k));
            hb(i, j, k) += _c_hb_psi_ea(ii, j, k) * _hb_psi_ea(ii, j, k);
          }
        }
      }
      return;
    }
    case Axis::XYZ::Y: {
      for (std::size_t i{0}; i < _nb; ++i) {
        for (std::size_t j{hNodeStartIndexMainAxis()};
             j < hNodeStartIndexMainAxis() + n(); ++j) {
          for (std::size_t k{0}; k < _na + 1; ++k) {
            auto jj{j - hNodeStartIndexMainAxis()};
            _ha_psi_eb(i, jj, k) =
                _coeff_b_h(jj) * _ha_psi_eb(i, jj, k) +
                _coeff_a_h(jj) * (eb(i, j + 1, k) - eb(i, j, k));
            ha(i, j, k) += _c_ha_psi_eb(i, jj, k) * _ha_psi_eb(i, jj, k);
          }
        }
      }

      for (std::size_t i{0}; i < _nb + 1; ++i) {
        for (std::size_t j{hNodeStartIndexMainAxis()};
             j < hNodeStartIndexMainAxis() + n(); ++j) {
          for (std::size_t k{0}; k < _na; ++k) {
            auto jj{j - hNodeStartIndexMainAxis()};
            _hb_psi_ea(i, jj, k) =
                _coeff_b_h(jj) * _hb_psi_ea(i, jj, k) +
                _coeff_a_h(jj) * (ea(i, j + 1, k) - ea(i, j, k));
            hb(i, j, k) += _c_hb_psi_ea(i, jj, k) * _hb_psi_ea(i, jj, k);
          }
        }
      }
      return;
    }
    case Axis::XYZ::Z: {
      for (std::size_t i{0}; i < _na + 1; ++i) {
        for (std::size_t j{0}; j < _nb; ++j) {
          for (std::size_t k{hNodeStartIndexMainAxis()};
               k < hNodeStartIndexMainAxis() + n(); ++k) {
            auto kk{k - hNodeStartIndexMainAxis()};
            _ha_psi_eb(i, j, kk) =
                _coeff_b_h(kk) * _ha_psi_eb(i, j, kk) +
                _coeff_a_h(i, j, kk) * (eb(i, j, k + 1) - eb(i, j, k));
            ha(i, j, k) += _c_ha_psi_eb(i, j, kk) * _ha_psi_eb(i, j, kk);
          }
        }
      }
      for (std::size_t i{0}; i < _na; ++i) {
        for (std::size_t j{0}; j < _nb + 1; ++j) {
          for (std::size_t k{hNodeStartIndexMainAxis()};
               k < hNodeStartIndexMainAxis() + n(); ++k) {
            auto kk{k - hNodeStartIndexMainAxis()};
            _hb_psi_ea(i, j, kk) =
                _coeff_b_h(kk) * _hb_psi_ea(i, j, kk) +
                _coeff_a_h(i, j, kk) * (ea(i, j, k + 1) - ea(i, j, k));
            hb(i, j, k) += _c_hb_psi_ea(i, j, kk) * _hb_psi_ea(i, j, kk);
          }
        }
      }
      return;
    }
    default:
      throw XFDTDPMLException("Invalid direction");
  }
}

void PML::correctCoefficientX() {
  _c_ea_psi_hb = xt::zeros<double>({n(), _na, _nb + 1});
  _ea_psi_hb = xt::zeros<double>({n(), _na, _nb + 1});

  auto& ceahb{calculationParamPtr()->fdtdCoefficient()->ceyhz()};

  _c_eb_psi_ha = xt::zeros<double>({n(), _na + 1, _nb});
  _eb_psi_ha = xt::zeros<double>({n(), _na + 1, _nb});

  auto& cebha{calculationParamPtr()->fdtdCoefficient()->cezhy()};

  for (std::size_t i{eNodeStartIndexMainAxis()};
       i < eNodeStartIndexMainAxis() + n(); ++i) {
    auto ii{i - eNodeStartIndexMainAxis()};
    for (std::size_t j{0}; j < _na; ++j) {
      for (std::size_t k{0}; k < _nb + 1; ++k) {
        _c_ea_psi_hb(ii, j, k) = ceahb(i, j, k) * _h_size(ii);
        ceahb(i, j, k) = ceahb(i, j, k) / _kappa_e(ii);
      }
    }

    for (std::size_t j{0}; j < _na + 1; ++j) {
      for (std::size_t k{0}; k < _nb; ++k) {
        _c_eb_psi_ha(ii, j, k) = cebha(i, j, k) * _h_size(ii);
        cebha(i, j, k) = cebha(i, j, k) / _kappa_e(ii);
      }
    }
  }

  _c_ha_psi_eb = xt::zeros<double>({n(), _na + 1, _nb});
  _ha_psi_eb = xt::zeros<double>({n(), _na + 1, _nb});

  auto& chaeb{calculationParamPtr()->fdtdCoefficient()->chyez()};

  _c_hb_psi_ea = xt::zeros<double>({n(), _na, _nb + 1});
  _hb_psi_ea = xt::zeros<double>({n(), _na, _nb + 1});

  auto& chbea{calculationParamPtr()->fdtdCoefficient()->chzey()};

  for (std::size_t i{hNodeStartIndexMainAxis()};
       i < hNodeStartIndexMainAxis() + n(); ++i) {
    auto ii{i - hNodeStartIndexMainAxis()};
    for (std::size_t j{0}; j < _na + 1; ++j) {
      for (std::size_t k{0}; k < _nb; ++k) {
        _c_ha_psi_eb(ii, j, k) = chaeb(i, j, k) * _e_size(ii);
        chaeb(i, j, k) = chaeb(i, j, k) / _kappa_h(ii);
      }
    }

    for (std::size_t j{0}; j < _na; ++j) {
      for (std::size_t k{0}; k < _nb + 1; ++k) {
        _c_hb_psi_ea(ii, j, k) = chbea(i, j, k) * _e_size(ii);
        chbea(i, j, k) = chbea(i, j, k) / _kappa_h(ii);
      }
    }
  }
}

void PML::correctCoefficientY() {
  _c_ea_psi_hb = xt::zeros<double>({_nb + 1, n(), _na});
  _ea_psi_hb = xt::zeros<double>({_nb + 1, n(), _na});
  auto& ceahb{calculationParamPtr()->fdtdCoefficient()->cezhx()};

  _c_eb_psi_ha = xt::zeros<double>({_nb, n(), _na + 1});
  _eb_psi_ha = xt::zeros<double>({_nb, n(), _na + 1});
  auto& cebha{calculationParamPtr()->fdtdCoefficient()->cexhz()};

  for (std::size_t i{0}; i < _nb + 1; ++i) {
    for (std::size_t j{eNodeStartIndexMainAxis()};
         j < eNodeStartIndexMainAxis() + n(); ++j) {
      auto jj{j - eNodeStartIndexMainAxis()};
      for (std::size_t k{0}; k < _na; ++k) {
        _c_ea_psi_hb(i, jj, k) = ceahb(i, j, k) * _h_size(jj);
        ceahb(i, j, k) = ceahb(i, j, k) / _kappa_e(jj);
      }
    }
  }

  for (std::size_t i{0}; i < _nb; ++i) {
    for (std::size_t j{eNodeStartIndexMainAxis()};
         j < eNodeStartIndexMainAxis() + n(); ++j) {
      auto jj{j - eNodeStartIndexMainAxis()};
      for (std::size_t k{0}; k < _na + 1; ++k) {
        _c_eb_psi_ha(i, jj, k) = cebha(i, j, k) * _h_size(jj);
        cebha(i, j, k) = cebha(i, j, k) / _kappa_e(jj);
      }
    }
  }

  _c_ha_psi_eb = xt::zeros<double>({_nb, n(), _na + 1});
  _ha_psi_eb = xt::zeros<double>({_nb, n(), _na + 1});
  auto& chaeb{calculationParamPtr()->fdtdCoefficient()->chzex()};

  _c_hb_psi_ea = xt::zeros<double>({_nb + 1, n(), _na});
  _hb_psi_ea = xt::zeros<double>({_nb + 1, n(), _na});
  auto& chbea{calculationParamPtr()->fdtdCoefficient()->chxez()};

  for (std::size_t i{0}; i < _nb; ++i) {
    for (std::size_t j{hNodeStartIndexMainAxis()};
         j < hNodeStartIndexMainAxis() + n(); ++j) {
      auto jj{j - hNodeStartIndexMainAxis()};
      for (std::size_t k{0}; k < _na + 1; ++k) {
        _c_ha_psi_eb(i, jj, k) = chaeb(i, j, k) * _e_size(jj);
        chaeb(i, j, k) = chaeb(i, j, k) / _kappa_h(jj);
      }
    }
  }

  for (std::size_t i{0}; i < _nb + 1; ++i) {
    for (std::size_t j{hNodeStartIndexMainAxis()};
         j < hNodeStartIndexMainAxis() + n(); ++j) {
      auto jj{j - hNodeStartIndexMainAxis()};
      for (std::size_t k{0}; k < _na; ++k) {
        _c_hb_psi_ea(i, jj, k) = chbea(i, j, k) * _e_size(jj);
        chbea(i, j, k) = chbea(i, j, k) / _kappa_h(jj);
      }
    }
  }
}

void PML::correctCoefficientZ() {
  _c_ea_psi_hb = xt::zeros<double>({_na, _nb + 1, n()});
  _ea_psi_hb = xt::zeros<double>({_na, _nb + 1, n()});
  auto& ceahb{calculationParamPtr()->fdtdCoefficient()->cexhy()};

  _c_eb_psi_ha = xt::zeros<double>({_na + 1, _nb, n()});
  _eb_psi_ha = xt::zeros<double>({_na + 1, _nb, n()});
  auto& cebha{calculationParamPtr()->fdtdCoefficient()->ceyhx()};

  for (std::size_t i{0}; i < _na; ++i) {
    for (std::size_t j{0}; j < _nb + 1; ++j) {
      for (std::size_t k{eNodeStartIndexMainAxis()};
           k < eNodeStartIndexMainAxis() + n(); ++k) {
        auto kk{k - eNodeStartIndexMainAxis()};
        _c_ea_psi_hb(i, j, kk) = ceahb(i, j, k) * _h_size(i, j, kk);
        ceahb(i, j, k) = ceahb(i, j, k) / _kappa_e(kk);
      }
    }
  }

  for (std::size_t i{0}; i < _na + 1; ++i) {
    for (std::size_t j{0}; j < _nb; ++j) {
      for (std::size_t k{eNodeStartIndexMainAxis()};
           k < eNodeStartIndexMainAxis() + n(); ++k) {
        auto kk{k - eNodeStartIndexMainAxis()};
        _c_eb_psi_ha(i, j, kk) = cebha(i, j, k) * _h_size(i, j, kk);
        cebha(i, j, k) = cebha(i, j, k) / _kappa_e(kk);
      }
    }
  }

  _c_ha_psi_eb = xt::zeros<double>({_na + 1, _nb, n()});
  _ha_psi_eb = xt::zeros<double>({_na + 1, _nb, n()});

  auto& chaeb{calculationParamPtr()->fdtdCoefficient()->chxey()};

  _c_hb_psi_ea = xt::zeros<double>({_na, _nb + 1, n()});
  _hb_psi_ea = xt::zeros<double>({_na, _nb + 1, n()});

  auto& chbea{calculationParamPtr()->fdtdCoefficient()->chyex()};

  for (std::size_t i{0}; i < _na + 1; ++i) {
    for (std::size_t j{0}; j < _nb; ++j) {
      for (std::size_t k{hNodeStartIndexMainAxis()};
           k < hNodeStartIndexMainAxis() + n(); ++k) {
        auto kk{k - hNodeStartIndexMainAxis()};
        _c_ha_psi_eb(i, j, kk) = chaeb(i, j, k) * _h_size(kk);
        chaeb(i, j, k) = chaeb(i, j, k) / _kappa_h(kk);
      }
    }
  }

  for (std::size_t i{0}; i < _na; ++i) {
    for (std::size_t j{0}; j < _nb+1; ++j) {
      for (std::size_t k{hNodeStartIndexMainAxis()};
           k < hNodeStartIndexMainAxis() + n(); ++k) {
        auto kk{k - hNodeStartIndexMainAxis()};
        _c_hb_psi_ea(i, j, kk) = chbea(i, j, k) * _h_size(kk);
        chbea(i, j, k) = chbea(i, j, k) / _kappa_h(kk);
      }
    }
  }
}

void PML::calRecursiveConvolutionCoeff() {
  double min_h_size{xt::amin(hSize())()};
  double sigma_max_e = calculateSigmaMax(min_h_size);
  auto rho_e{calculateRhoE(n(), hSize())};
  auto sigma_e{calculateSigma(sigma_max_e, rho_e, _order)};
  auto alpha_e{calculateAlpha(_alpha_min, _alpha_max, rho_e)};

  _kappa_e = calculateKappa(_kappa_max, rho_e, _order);

  _coeff_b_e = calculateCoefficientB(sigma_e, _kappa_e, alpha_e,
                                     calculationParamPtr()->timeParam()->dt(),
                                     constant::EPSILON_0);

  _coeff_a_e =
      calculateCoefficientA(_coeff_b_e, sigma_e, _kappa_e, alpha_e, hSize());

  double e2m{constant::MU_0 / constant::EPSILON_0};
  auto min_e_size{xt::amin(eSize())()};
  double sigma_max_m{calculateSigmaMax(min_e_size) * e2m};
  auto rho_m{calculateRhoM(n(), eSize())};
  auto sigma_m{calculateSigma(sigma_max_m, rho_m, _order)};
  auto alpha_m{calculateAlpha(_alpha_min, _alpha_max, rho_m) * e2m};

  _kappa_h = {calculateKappa(_kappa_max, rho_m, _order)};

  _coeff_b_h = calculateCoefficientB(sigma_m, _kappa_h, alpha_m,
                                     calculationParamPtr()->timeParam()->dt(),
                                     constant::MU_0);

  _coeff_a_h =
      calculateCoefficientA(_coeff_b_h, sigma_m, _kappa_h, alpha_m, eSize());
}

double PML::calculateSigmaMax(double dl) const {
  return _sigma_ratio * (_order + 1) / (150 * constant::PI * dl);
}

xt::xarray<double> PML::calculateRhoE(std::size_t n,
                                      const xt::xarray<double>& size) const {
  assert(n == size.size());
  auto interval{size - 0.75 * size};
  xt::xarray<double> d;
  d.resize({n});
  d(0) = interval(0);
  double sum{0};
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

xt::xarray<double> PML::calculateRhoM(std::size_t n,
                                      const xt::xarray<double>& size) const {
  assert(n == size.size());
  auto interval{size - 0.25 * size};
  xt::xarray<double> d;
  d.resize({n});
  d(0) = interval(0);
  double sum{0};
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

xt::xarray<double> PML::calculateSigma(double sigma_max,
                                       const xt::xarray<double>& rho,
                                       std::size_t order) const {
  return sigma_max * xt::pow(rho, _order);
}

xt::xarray<double> PML::calculateKappa(double kappa_max,
                                       const xt::xarray<double>& rho,
                                       std::size_t order) const {
  return 1 + (kappa_max - 1) * xt::pow(rho, _order);
}

xt::xarray<double> PML::calculateAlpha(double alpha_min, double alpha_max,
                                       const xt::xarray<double>& rho) const {
  return alpha_min + (alpha_max - alpha_min) * (1 - rho);
}

xt::xarray<double> PML::calculateCoefficientA(
    const xt::xarray<double>& b, const xt::xarray<double>& sigma,
    const xt::xarray<double>& kappa, const xt::xarray<double>& alpha,
    const xt::xarray<double>& dl) const {
  return (b - 1) * sigma / ((sigma + kappa * alpha) * kappa * dl);
}

xt::xarray<double> PML::calculateCoefficientB(const xt::xarray<double>& sigma,
                                              const xt::xarray<double>& kappa,
                                              const xt::xarray<double>& alpha,
                                              double dt,
                                              double constant) const {
  return xt::exp((-dt / constant) * (sigma / kappa + alpha));
}

xt::xarray<double> PML::calculateCoeffPsi(const xt::xarray<double>& coeff,
                                          const xt::xarray<double>& kappa,
                                          const xt::xarray<double>& dl) const {
  return coeff * dl;
}

}  // namespace xfdtd
