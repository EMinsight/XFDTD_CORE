#ifndef _XFDTD_LIB_FLUX_H_
#define _XFDTD_LIB_FLUX_H_

#include <cmath>
#include <complex>
#include <cstdlib>
#include <memory>
#include <tuple>

#include "xfdtd/coordinate_system/coordinate_system.h"
#include "xfdtd/electromagnetic_field/electromagnetic_field.h"
#include "xfdtd/exception/exception.h"
#include "xfdtd/monitor/monitor.h"
#include "xfdtd/shape/shape.h"
#include "xfdtd/util/constant.h"

namespace xfdtd {

class Flux : public Monitor {
 public:
  Flux(std::string name, std::unique_ptr<Shape> shape,
       Axis::Direction direction);

  Flux(const Flux&) = delete;

  Flux(Flux&&) noexcept = default;

  Flux& operator=(const Flux&) = delete;

  Flux& operator=(Flux&&) noexcept = default;

  ~Flux() override = default;

  Axis::Direction direction() const;

  Axis::XYZ mainAxis() const;

 private:
  Axis::Direction _direction;
  Axis::XYZ _main_axis;
};

Flux::Flux(std::string name, std::unique_ptr<Shape> shape,
           Axis::Direction direction)
    : Monitor(std::move(shape), std::move(name)), _direction{direction} {}

class FrequencyDomainFlux : public Flux {
 public:
  FrequencyDomainFlux(std::string name, std::unique_ptr<Shape> shape,
                      Axis::Direction direction,
                      xt::xarray<double> frequencies);

  FrequencyDomainFlux(const FrequencyDomainFlux&) = delete;

  FrequencyDomainFlux(FrequencyDomainFlux&&) noexcept = default;

  FrequencyDomainFlux& operator=(const FrequencyDomainFlux&) = delete;

  FrequencyDomainFlux& operator=(FrequencyDomainFlux&&) noexcept = default;

  ~FrequencyDomainFlux() override = default;

  void init(std::shared_ptr<const GridSpace> grid_space,
            std::shared_ptr<const CalculationParam> calculation_param,
            std::shared_ptr<const EMF> emf) override;

  void update() override;

  std::tuple<xt::xarray<std::complex<double>>, xt::xarray<std::complex<double>>,
             xt::xarray<std::complex<double>>, xt::xarray<std::complex<double>>>
  vectorPotential(const Vector& origin, double theta, double phi) const;

 private:
  xt::xarray<double> _frequencies;

  xt::xarray<std::complex<double>> _ja, _jb, _ma, _mb;
  EMF::Field _ea, _eb, _ha, _hb;
  xt::xarray<std::complex<double>> _transform_e, _transform_h;

  void initDFT();

  double waveNumber(std::size_t f) const;

  double delay(const Vector& origin, std::size_t i, std::size_t j,
               std::size_t k, double sin_t_cos_p, double sin_t_sin_p,
               double cos_t) const;

  double area(std::size_t i, std::size_t j, std::size_t k) const;
};

FrequencyDomainFlux::FrequencyDomainFlux(std::string name,
                                         std::unique_ptr<Shape> shape,
                                         Axis::Direction direction,
                                         xt::xarray<double> frequencies)
    : Flux{std::move(name), std::move(shape), direction},
      _frequencies{std::move(frequencies)} {}

void FrequencyDomainFlux::init(
    std::shared_ptr<const GridSpace> grid_space,
    std::shared_ptr<const CalculationParam> calculation_param,
    std::shared_ptr<const EMF> emf) {
  defaultInit(grid_space, calculation_param, emf);

  auto box{gridBoxPtr()};
  auto is{box->origin().i()};
  auto js{box->origin().j()};
  auto ks{box->origin().k()};
  auto ie{box->end().i()};
  auto je{box->end().j()};
  auto ke{box->end().k()};

  _ja = xt::zeros<std::complex<double>>(
      {_frequencies.size(), ie - is, je - js, ke - ks});
  _jb = xt::zeros<std::complex<double>>(
      {_frequencies.size(), ie - is, je - js, ke - ks});

  _ma = xt::zeros<std::complex<double>>(
      {_frequencies.size(), ie - is, je - js, ke - ks});
  _mb = xt::zeros<std::complex<double>>(
      {_frequencies.size(), ie - is, je - js, ke - ks});

  switch (direction()) {
    case Axis::Direction::XP:
      _ea = EMF::Field::EY;
      _eb = EMF::Field::EZ;
      _ha = EMF::Field::HY;
      _hb = EMF::Field::HZ;
      break;
    case Axis::Direction::YP:
      _ea = EMF::Field::EZ;
      _eb = EMF::Field::EX;
      _ha = EMF::Field::HZ;
      _hb = EMF::Field::HX;
      break;
    case Axis::Direction::ZP:
      _ea = EMF::Field::EX;
      _eb = EMF::Field::EY;
      _ha = EMF::Field::HX;
      _hb = EMF::Field::HY;
      break;
    case Axis::Direction::XN:
      _ea = EMF::Field::EZ;
      _eb = EMF::Field::EY;
      _ha = EMF::Field::HZ;
      _hb = EMF::Field::HY;
      break;
    case Axis::Direction::YN:
      _ea = EMF::Field::EX;
      _eb = EMF::Field::EZ;
      _ha = EMF::Field::HX;
      _hb = EMF::Field::HZ;
      break;
    case Axis::Direction::ZN:
      _ea = EMF::Field::EY;
      _eb = EMF::Field::EX;
      _ha = EMF::Field::HY;
      _hb = EMF::Field::HX;
      break;
    default:
      throw XFDTDException("Invalid direction");
  }

  initDFT();
}

void FrequencyDomainFlux::update() {
  auto box{gridBoxPtr()};
  const auto is{box->origin().i()};
  const auto js{box->origin().j()};
  const auto ks{box->origin().k()};
  const auto ie{box->end().i()};
  const auto je{box->end().j()};
  const auto ke{box->end().k()};

  auto emf{emfPtr()};
  auto t{calculationParamPtr()->timeParam()->currentTimeStep()};

  for (auto i{is}; i < ie; ++i) {
    for (auto j{js}; j < je; ++j) {
      for (auto k{ks}; k < ke; ++k) {
        auto ma{emf->fieldFaceCenter(i, j, k, _ea, mainAxis())};
        auto mb{-1 * emf->fieldFaceCenter(i, j, k, _eb, mainAxis())};
        auto ja{-1 * emf->fieldFaceCenter(i, j, k, _ha, mainAxis())};
        auto jb{emf->fieldFaceCenter(i, j, k, _hb, mainAxis())};

        for (auto f{0}; f < _frequencies.size(); ++f) {
          _ma(f, i, j, k) += ma * _transform_e(f, t);
          _mb(f, i, j, k) += mb * _transform_e(f, t);
          _ja(f, i, j, k) += ja * _transform_h(f, t);
          _jb(f, i, j, k) += jb * _transform_h(f, t);
        }
      }
    }
  }
}

std::tuple<xt::xarray<std::complex<double>>, xt::xarray<std::complex<double>>,
           xt::xarray<std::complex<double>>, xt::xarray<std::complex<double>>>
FrequencyDomainFlux::vectorPotential(const Vector& origin, double theta,
                                     double phi) const {
  using namespace std::complex_literals;

  xt::xarray<std::complex<double>> a_a{
      xt::zeros<std::complex<double>>({_frequencies.size()})};
  xt::xarray<std::complex<double>> a_b{
      xt::zeros<std::complex<double>>({_frequencies.size()})};
  xt::xarray<std::complex<double>> f_a{
      xt::zeros<std::complex<double>>({_frequencies.size()})};
  xt::xarray<std::complex<double>> f_b{
      xt::zeros<std::complex<double>>({_frequencies.size()})};

  auto sin_t_cos_p{std::sin(theta) * std::cos(phi)};
  auto sin_t_sin_p{std::sin(theta) * std::sin(phi)};
  auto cos_t{std::cos(theta)};

  auto box{gridBoxPtr()};
  const auto is{box->origin().i()};
  const auto js{box->origin().j()};
  const auto ks{box->origin().k()};
  const auto ie{box->end().i()};
  const auto je{box->end().j()};
  const auto ke{box->end().k()};
  for (auto f{0}; f < _frequencies.size(); ++f) {
    auto wave_number{waveNumber(f)};
    for (auto i{is}; i < ie; ++i) {
      for (auto j{js}; j < je; ++j) {
        for (auto k{ks}; k < ke; ++k) {
          auto ds{area(i, j, k)};
          auto phase_shift{std::exp(
              1i * wave_number *
              delay(origin, i, j, k, sin_t_cos_p, sin_t_sin_p, cos_t))};
          a_a(f) += ds * phase_shift * _ja(f, i, j, k);
          a_b(f) += ds * phase_shift * _jb(f, i, j, k);
          f_a(f) += ds * phase_shift * _ma(f, i, j, k);
          f_b(f) += ds * phase_shift * _mb(f, i, j, k);
        }
      }
    }
  }

  return {a_a, a_b, f_a, f_b};
}

void FrequencyDomainFlux::initDFT() {
  using namespace std::complex_literals;
  const auto total_time_step{
      calculationParamPtr()->timeParam()->endTimeStep() -
      calculationParamPtr()->timeParam()->startTimeStep() + 1};
  const auto dt{calculationParamPtr()->timeParam()->dt()};
  _transform_e =
      xt::zeros<std::complex<double>>({_frequencies.size(), total_time_step});
  _transform_h =
      xt::zeros<std::complex<double>>({_frequencies.size(), total_time_step});
  for (auto f{0}; f < _frequencies.size(); ++f) {
    for (std::size_t t{0}; t < total_time_step; ++t) {
      _transform_e(f, t) =
          dt * std::exp(-1i * 2.0 * constant::PI * (_frequencies(f) * t * dt));
      _transform_h(f, t) = dt * std::exp(-1i * 2.0 * constant::PI *
                                         (_frequencies(f) * (t - 0.5) * dt));
    }
  }
}

double FrequencyDomainFlux::waveNumber(std::size_t f) const {
  return 2 * constant::PI * _frequencies(f) / constant::C_0;
}

double FrequencyDomainFlux::delay(const Vector& origin, std::size_t i,
                                  std::size_t j, std::size_t k,
                                  double sin_t_cos_p, double sin_t_sin_p,
                                  double cos_t) const {
  auto g{gridSpacePtr()};
  auto main_axis{Axis::formDirectionToXYZ(direction())};
  const auto& h_node_x{g->hNodeX()};
  const auto& h_node_y{g->hNodeY()};
  const auto& h_node_z{g->hNodeZ()};
  const auto& e_node_x{g->eNodeX()};
  const auto& e_node_y{g->eNodeY()};
  const auto& e_node_z{g->eNodeZ()};
  double res{0};
  switch (main_axis) {
    case Axis::XYZ::X:
      res = std::abs(origin.x() - h_node_x(i)) * sin_t_cos_p +
            std::abs(origin.y() - e_node_y(j)) * sin_t_sin_p +
            std::abs(origin.z() - e_node_z(k)) * cos_t;
      break;
    case Axis::XYZ::Y:
      res = std::abs(origin.x() - e_node_x(i)) * sin_t_cos_p +
            std::abs(origin.y() - h_node_y(j)) * sin_t_sin_p +
            std::abs(origin.z() - e_node_z(k)) * cos_t;
      break;
    case Axis::XYZ::Z:
      res = std::abs(origin.x() - e_node_x(i)) * sin_t_cos_p +
            std::abs(origin.y() - e_node_y(j)) * sin_t_sin_p +
            std::abs(origin.z() - h_node_z(k)) * cos_t;
      break;
    default:
      throw XFDTDException("Invalid main axis");
  }
  return res;
}

double FrequencyDomainFlux::area(std::size_t i, std::size_t j,
                                 std::size_t k) const {
  auto g{gridSpacePtr()};
  auto main_axis{Axis::formDirectionToXYZ(direction())};
  const auto& h_size_x{g->hSizeX()};
  const auto& h_size_y{g->hSizeY()};
  const auto& h_size_z{g->hSizeZ()};
  double res{0};
  switch (main_axis) {
    case Axis::XYZ::X:
      res = h_size_y(j) * h_size_z(k);
      break;
    case Axis::XYZ::Y:
      res = h_size_x(i) * h_size_z(k);
      break;
    case Axis::XYZ::Z:
      res = h_size_x(i) * h_size_y(j);
      break;
    default:
      throw XFDTDException("Invalid main axis");
  }
  return res;
}

}  // namespace xfdtd

#endif  // _XFDTD_LIB_FLUX_H_
