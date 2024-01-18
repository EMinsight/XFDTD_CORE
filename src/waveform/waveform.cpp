#include <xfdtd/util/constant.h>
#include <xfdtd/waveform/waveform.h>

#include <utility>
#include <xtensor.hpp>

namespace xfdtd {

Waveform::Waveform(std::function<double(double)> func, double amplitude)
    : _func{std::move(func)}, _amplitude{amplitude} {}

double Waveform::operator()(double t) { return _func(t); }

std::function<double(double)> Waveform::func() const { return _func; }

double Waveform::amplitude() const { return _amplitude; }

const xt::xarray<double>& Waveform::value() const { return _value; }

const xt::xarray<double>& Waveform::time() const { return _time; }

xt::xarray<double>& Waveform::value() { return _value; }

xt::xarray<double>& Waveform::time() { return _time; }

void Waveform::init(xt::xarray<double> time) {
  _value = xt::zeros<double>(time.shape());
  for (auto i = 0; i < time.size(); ++i) {
    _value(i) = amplitude() * _func(time(i));
  }
  _time = std::move(time);
}

void Waveform::setAmplitude(double amplitude) { _amplitude = amplitude; }

Waveform Waveform::sine(double frequency, double amplitude) {
  return Waveform{[frequency](double t) {
                    return std::sin(2 * constant::PI * frequency * t);
                  },
                  amplitude};
}

Waveform Waveform::cosine(double frequency, double amplitude) {
  return Waveform{[frequency](double t) {
                    return std::cos(2 * constant::PI * frequency * t);
                  },
                  amplitude};
}

Waveform Waveform::square(double frequency, double amplitude) {
  return Waveform(
      [frequency](double t) {
        return std::copysign(1.0, std::sin(2 * constant::PI * frequency * t));
      },
      amplitude);
}

Waveform Waveform::triangle(double frequency, double amplitude) {
  return Waveform([frequency](double t) {
    return std::asin(std::sin(2 * constant::PI * frequency * t)) * 2 /
           constant::PI;
  });
}

Waveform Waveform::sawtooth(double frequency, double amplitude) {
  return Waveform(
      [frequency](double t) { return (2 * std::fmod(frequency * t, 1.0) - 1); },
      amplitude);
}

Waveform Waveform::gaussian(double tau, double t_0, double amplitude) {
  return Waveform(
      [tau, t_0](double t) { return std::exp(-std::pow((t - t_0) / tau, 2)); },
      amplitude);
}

Waveform Waveform::cosineModulatedGaussian(double tau, double t_0,
                                           double frequency, double amplitude) {
  return Waveform(
      [tau, t_0, frequency](double t) {
        return std::exp(-std::pow((t - t_0) / tau, 2)) *
               std::cos(2 * constant::PI * frequency * (t - t_0));
      },
      amplitude);
}

Waveform Waveform::step(double t_0, double amplitude) {
  return Waveform([t_0](double t) { return static_cast<double>(t >= t_0); },
                  amplitude);
}

}  // namespace xfdtd
