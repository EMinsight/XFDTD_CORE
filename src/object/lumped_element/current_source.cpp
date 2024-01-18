#include <xfdtd/object/lumped_element/current_source.h>

#include <xtensor.hpp>

namespace xfdtd {

CurrentSource::CurrentSource(std::string name, std::unique_ptr<Cube> cube,
                             Axis::Direction direction, double resistance,
                             std::unique_ptr<Waveform> waveform,
                             std::unique_ptr<Material> material)
    : LumpedElement{std::move(name), std::move(cube),
                    Axis::formDirectionToXYZ(direction), std::move(material)},
      _direction{direction},
      _resistance{resistance},
      _waveform{std::move(waveform)} {
  if (_resistance == 0) {
    _resistance = 1e-20;
  }
}

std::string CurrentSource::toString() const {
  // return "CurrentSource{name: " + name() + ", shape: " + shape()->toString()
  // +
  //        ", material: " + material()->toString() +
  //        ", direction: " + Axis::toString(_direction) +
  //        ", resistance: " + std::to_string(_resistance) +
  //        ", waveform: " + _waveform->toString() + "}";
  return "";
}

Axis::Direction CurrentSource::direction() const { return _direction; }

double CurrentSource::resistance() const { return _resistance; }

const std::unique_ptr<Waveform> &CurrentSource::waveform() const {
  return _waveform;
}

void CurrentSource::init(std::shared_ptr<const GridSpace> grid_space,
                         std::shared_ptr<CalculationParam> calculation_param,
                         std::shared_ptr<EMF> emf) {
  LumpedElement::init(std::move(grid_space), std::move(calculation_param),
                      std::move(emf));

  auto rf{[](double r, std::size_t na, std::size_t nb, std::size_t nc) {
    return r * na * nb / nc;
  }};

  auto cf{[](double c, std::size_t na, std::size_t nb, std::size_t nc) {
    return c / (na * nb);
  }};

  auto dx_dy_dz{[](const xt::xarray<double> &x, const xt::xarray<double> &y,
                   const xt::xarray<double> &z, auto &&x_range, auto &&y_range,
                   auto &&z_range) {
    return xt::meshgrid(xt::view(x, x_range), xt::view(y, y_range),
                        xt::view(z, z_range));
  }};

  _resistance_factor = rf(_resistance, nodeCountSubAxisA(), nodeCountSubAxisB(),
                          nodeCountMainAxis());
  _current_amplitude_factor = cf(_waveform->amplitude(), nodeCountSubAxisA(),
                                 nodeCountSubAxisB(), nodeCountMainAxis());

  auto g{gridSpacePtr()};

  if (xyz() == Axis::XYZ::X) {
    if (_direction == Axis::Direction::XN) {
      _current_amplitude_factor *= -1;
    }

    auto [dx, dy, dz] = dx_dy_dz(g->eSizeX(), g->hSizeY(), g->hSizeZ(),
                                 rangeX(), rangeY(), rangeZ());

    _da = dy;
    _db = dz;
    _dc = dx;
  }

  if (xyz() == Axis::XYZ::Y) {
    if (_direction == Axis::Direction::YN) {
      _current_amplitude_factor *= -1;
    }

    auto [dx, dy, dz] = dx_dy_dz(g->hSizeX(), g->eSizeY(), g->hSizeZ(),
                                 rangeX(), rangeY(), rangeZ());

    _da = dz;
    _db = dx;
    _dc = dy;
  }

  if (xyz() == Axis::XYZ::Z) {
    if (_direction == Axis::Direction::ZN) {
      _current_amplitude_factor *= -1;
    }

    auto [dx, dy, dz] = dx_dy_dz(g->hSizeX(), g->hSizeY(), g->eSizeZ(),
                                 rangeX(), rangeY(), rangeZ());

    _da = dx;
    _db = dy;
    _dc = dz;
  }

  auto dt{calculationParamPtr()->timeParam()->dt()};
  _beta = (dt * _dc) / (_resistance_factor * _da * _db);
  _waveform->init(calculationParamPtr()->timeParam()->hTime());
}

void CurrentSource::correctUpdateCoefficient() {
  auto dt{calculationParamPtr()->timeParam()->dt()};
  auto func{[this, dt](xt::xarray<double> &cece, xt::xarray<double> &cecha,
                       xt::xarray<double> &cechb, const xt::xarray<double> &eps,
                       const xt::xarray<double> &sigma) {
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
    _coff_i = -2 * dt /
              ((2 * eps_view + dt * sigma_view + _beta) * (_db * _da)) *
              _current_amplitude_factor;
  }};

  auto calc_param{calculationParamPtr()};

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

void CurrentSource::correctE() {
  auto func{[this](xt::xarray<double> &e) {
    auto range_x{rangeX()};
    auto range_y{rangeY()};
    auto range_z{rangeZ()};
    auto e_view{xt::view(e, range_x, range_y, range_z)};
    e_view +=
        _coff_i *
        _waveform
            ->value()[calculationParamPtr()->timeParam()->currentTimeStep()];
  }};

  auto emf{emfPtr()};
  if (xyz() == Axis::XYZ::X) {
    func(emf->ex());
    return;
  }

  if (xyz() == Axis::XYZ::Y) {
    func(emf->ey());
    return;
  }

  if (xyz() == Axis::XYZ::Z) {
    func(emf->ez());
    return;
  }
}

void CurrentSource::correctH() {}

}  // namespace xfdtd
