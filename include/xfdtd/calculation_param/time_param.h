#ifndef __XFDTD_CORE_TIME_PARAM_H__
#define __XFDTD_CORE_TIME_PARAM_H__

#include <xfdtd/common/type_define.h>
#include <xfdtd/exception/exception.h>

namespace xfdtd {

class XFDTDTimeParamException : public XFDTDException {
 public:
  explicit XFDTDTimeParamException(
      std::string message = "XFDTD Time Parameter Exception")
      : XFDTDException(std::move(message)) {}
};

class TimeParam {
 public:
  static constexpr Real DEFAULT_CFL{0.98};

  static Real calculateDt(Real cfl, Real dx, Real dy, Real dz);

  static Real calculateDt(Real cfl, Real dx, Real dy);

  static Real calculateDt(Real cfl, Real dz);

  TimeParam(Real dt, std::size_t size, std::size_t start_time_step = 0);

  explicit TimeParam(Real cfl = DEFAULT_CFL);

  Real dt() const;

  Real cfl() const;

  std::size_t startTimeStep() const;

  std::size_t size() const;

  std::size_t endTimeStep() const;

  std::size_t currentTimeStep() const;

  std::size_t remainingTimeStep() const;

  virtual void nextStep();

  virtual void reset();

  void setDt(Real dt);

  void setTimeParamRunRange(std::size_t end_time_step,
                            std::size_t start_time_step = 0);

  Array1D<Real> eTime() const;

  Array1D<Real> hTime() const;

 private:
  Real _cfl;
  Real _dt{};
  std::size_t _start_time_step{};
  std::size_t _size{};
  std::size_t _current_time_step{};
};

}  // namespace xfdtd

#endif  // __XFDTD_CORE_TIME_PARAM_H__
