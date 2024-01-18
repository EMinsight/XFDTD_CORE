#ifndef _XFDTD_LIB_WAVEFORM_H_
#define _XFDTD_LIB_WAVEFORM_H_

#include <functional>
#include <xtensor/xarray.hpp>

namespace xfdtd {

class Waveform {
 public:
  explicit Waveform(std::function<double(double)> func, double amplitude = 1.0);

  Waveform(const Waveform&) = default;

  Waveform(Waveform&&) noexcept = default;

  Waveform& operator=(const Waveform&) = default;

  Waveform& operator=(Waveform&&) noexcept = default;

  ~Waveform() = default;

  static Waveform sine(double frequency, double amplitude = 1.0);

  static Waveform cosine(double frequency, double amplitude = 1.0);

  static Waveform square(double frequency, double amplitude = 1.0);

  static Waveform triangle(double frequency, double amplitude = 1.0);

  static Waveform sawtooth(double frequency, double amplitude = 1.0);

  static Waveform gaussian(double tau, double t_0, double amplitude = 1.0);

  static Waveform cosineModulatedGaussian(double tau, double t_0,
                                          double frequency,
                                          double amplitude = 1.0);

  static Waveform step(double t_0, double amplitude = 1.0);

  double operator()(double t);

  std::function<double(double)> func() const;

  double amplitude() const;

  const xt::xarray<double>& value() const;

  const xt::xarray<double>& time() const;

  xt::xarray<double>& value();

  xt::xarray<double>& time();

  void init(xt::xarray<double> time);

  void setAmplitude(double amplitude);

 private:
  std::function<double(double)> _func;

  double _amplitude;

  xt::xarray<double> _time;
  xt::xarray<double> _value;
};

}  // namespace xfdtd

#endif  // _XFDTD_LIB_WAVEFORM_H_
