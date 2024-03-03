#ifndef _XFDTD_LIB_TFSF_1D_H_
#define _XFDTD_LIB_TFSF_1D_H_

#include "xfdtd/waveform_source/tfsf.h"

namespace xfdtd {

class TFSF1D : public TFSF {
 public:
  TFSF1D(std::size_t z, bool forward, std::unique_ptr<Waveform> waveform);

  TFSF1D(TFSF1D &&) noexcept = default;

  TFSF1D &operator=(TFSF1D &&) noexcept = default;

  ~TFSF1D() override = default;

  void init(std::shared_ptr<GridSpace> grid_space,
            std::shared_ptr<CalculationParam> calculation_param,
            std::shared_ptr<EMF> emf) override;

  void correctE() override;

  void correctH() override;

  bool forward() const;

  std::size_t start() const;

  std::unique_ptr<Corrector> generateCorrector(
      const Divider::IndexTask &task) override;

 private:
  bool _forward;
  std::size_t _start;
};

}  // namespace xfdtd

#endif  // _XFDTD_LIB_TFSF_1D_H_
