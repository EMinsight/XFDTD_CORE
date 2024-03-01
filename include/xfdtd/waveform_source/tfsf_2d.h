#ifndef _XFDTD_LIB_TFSF_2D_H_
#define _XFDTD_LIB_TFSF_2D_H_

#include "xfdtd/waveform_source/tfsf.h"

namespace xfdtd {

class TFSF2D : public TFSF {
 public:
  TFSF2D(std::size_t distance_x, std::size_t distance_y, double phi,
         std::unique_ptr<Waveform> waveform);

  ~TFSF2D() override = default;

  void init(std::shared_ptr<GridSpace> grid_space,
            std::shared_ptr<CalculationParam> calculation_param,
            std::shared_ptr<EMF> emf) override;

  void correctE() override;

  void correctH() override;

  std::unique_ptr<Corrector> generateCorrector(
      const Divider::IndexTask &task) override;

 private:
};

}  // namespace xfdtd

#endif  // _XFDTD_LIB_TFSF_2D_H_
