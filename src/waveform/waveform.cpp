#include <xfdtd/common/constant.h>
#include <xfdtd/common/type_define.h>
#include <xfdtd/waveform/waveform.h>

#include <utility>
#include <xtensor.hpp>

namespace xfdtd {

Waveform::Waveform(std::function<Real(Real)> func, Real amplitude)
    : _func{std::move(func)}, _amplitude{amplitude} {}

Real Waveform::operator()(Real t) { return _func(t); }

std::function<Real(Real)> Waveform::func() const { return _func; }

Real Waveform::amplitude() const { return _amplitude; }

const Array1D<Real>& Waveform::value() const { return _value; }

const Array1D<Real>& Waveform::time() const { return _time; }

Array1D<Real>& Waveform::value() { return _value; }

Array1D<Real>& Waveform::time() { return _time; }

void Waveform::init(Array1D<Real> time) {
  _value = xt::zeros<Real>(time.shape());
  for (auto i = 0; i < time.size(); ++i) {
    _value(i) = amplitude() * _func(time(i));
  }
  _time = std::move(time);
}

void Waveform::setAmplitude(Real amplitude) { _amplitude = amplitude; }

Waveform Waveform::sine(Real frequency, Real amplitude) {
  return Waveform{[frequency](Real t) {
                    return std::sin(2 * constant::PI * frequency * t);
                  },
                  amplitude};
}

Waveform Waveform::cosine(Real frequency, Real amplitude) {
  return Waveform{[frequency](Real t) {
                    return std::cos(2 * constant::PI * frequency * t);
                  },
                  amplitude};
}

Waveform Waveform::square(Real frequency, Real amplitude) {
  return Waveform(
      [frequency](Real t) {
        return std::copysign(1.0, std::sin(2 * constant::PI * frequency * t));
      },
      amplitude);
}

Waveform Waveform::triangle(Real frequency, Real amplitude) {
  return Waveform([frequency](Real t) {
    return std::asin(std::sin(2 * constant::PI * frequency * t)) * 2 /
           constant::PI;
  });
}

Waveform Waveform::sawtooth(Real frequency, Real amplitude) {
  return Waveform(
      [frequency](Real t) { return (2 * std::fmod(frequency * t, 1.0) - 1); },
      amplitude);
}

Waveform Waveform::gaussian(Real tau, Real t_0, Real amplitude) {
  return Waveform(
      [tau, t_0](Real t) { return std::exp(-std::pow((t - t_0) / tau, 2)); },
      amplitude);
}

Waveform Waveform::cosineModulatedGaussian(Real tau, Real t_0,
                                           Real frequency, Real amplitude) {
  return Waveform(
      [tau, t_0, frequency](Real t) {
        return std::exp(-std::pow((t - t_0) / tau, 2)) *
               std::cos(2 * constant::PI * frequency * (t - t_0));
      },
      amplitude);
}

Waveform Waveform::step(Real t_0, Real amplitude) {
  return Waveform([t_0](Real t) { return static_cast<Real>(t >= t_0); },
                  amplitude);
}

}  // namespace xfdtd
