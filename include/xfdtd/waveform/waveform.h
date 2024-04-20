#ifndef _XFDTD_CORE_WAVEFORM_H_
#define _XFDTD_CORE_WAVEFORM_H_

#include <xfdtd/common/type_define.h>

#include <functional>
#include <memory>

namespace xfdtd {

class Waveform {
 public:
  explicit Waveform(std::function<Real(Real)> func, Real amplitude = 1.0);

  Waveform(const Waveform&) = default;

  Waveform(Waveform&&) noexcept = default;

  Waveform& operator=(const Waveform&) = default;

  Waveform& operator=(Waveform&&) noexcept = default;

  ~Waveform() = default;

  static auto sine(Real frequency, Real amplitude = 1.0)
      -> std::unique_ptr<Waveform>;

  static auto cosine(Real frequency, Real amplitude = 1.0)
      -> std::unique_ptr<Waveform>;

  static auto square(Real frequency, Real amplitude = 1.0)
      -> std::unique_ptr<Waveform>;

  static auto triangle(Real frequency, Real amplitude = 1.0)
      -> std::unique_ptr<Waveform>;

  static auto sawtooth(Real frequency, Real amplitude = 1.0)
      -> std::unique_ptr<Waveform>;

  static auto gaussian(Real tau, Real t_0, Real amplitude = 1.0)
      -> std::unique_ptr<Waveform>;

  static auto cosineModulatedGaussian(Real tau, Real t_0, Real frequency,
                                      Real amplitude = 1.0)
      -> std::unique_ptr<Waveform>;

  static auto step(Real t_0, Real amplitude = 1.0) -> std::unique_ptr<Waveform>;

  Real operator()(Real t);

  std::function<Real(Real)> func() const;

  Real amplitude() const;

  const Array1D<Real>& value() const;

  const Array1D<Real>& time() const;

  Array1D<Real>& value();

  Array1D<Real>& time();

  void init(Array1D<Real> time);

  void setAmplitude(Real amplitude);

 private:
  std::function<Real(Real)> _func;

  Real _amplitude;

  Array1D<Real> _time;
  Array1D<Real> _value;
};

}  // namespace xfdtd

#endif  // _XFDTD_CORE_WAVEFORM_H_
