#ifndef _XFDTD_CORE_TFSF_3D_H_
#define _XFDTD_CORE_TFSF_3D_H_

#include "xfdtd/waveform_source/tfsf.h"

namespace xfdtd {

class TFSF3D : public TFSF {
 public:
  TFSF3D(std::size_t x, std::size_t y, std::size_t z, Real theta, Real phi,
         Real psi, std::unique_ptr<Waveform> waveform);

  TFSF3D(TFSF3D &&) noexcept = default;

  TFSF3D &operator=(TFSF3D &&) noexcept = default;

  ~TFSF3D() override = default;

  void init(std::shared_ptr<GridSpace> grid_space,
            std::shared_ptr<CalculationParam> calculation_param,
            std::shared_ptr<EMF> emf) override;

  std::unique_ptr<Corrector> generateCorrector(
      const IndexTask &task) override;

 private:
};

}  // namespace xfdtd

#endif  // _XFDTD_CORE_TFSF_3D_H_
