#ifndef _XFDTD_CORE_TFSF_1D_H_
#define _XFDTD_CORE_TFSF_1D_H_

#include <xfdtd/common/type_define.h>
#include <xfdtd/waveform_source/tfsf.h>

namespace xfdtd {

class TFSF1D : public TFSF {
 public:
  TFSF1D(Index z, bool forward, std::unique_ptr<Waveform> waveform);

  TFSF1D(TFSF1D &&) noexcept = default;

  TFSF1D &operator=(TFSF1D &&) noexcept = default;

  ~TFSF1D() override = default;

  auto forward() const -> bool;

  auto start() const -> Index;

  auto init(std::shared_ptr<GridSpace> grid_space,
            std::shared_ptr<CalculationParam> calculation_param,
            std::shared_ptr<EMF> emf) -> void override;

  auto generateCorrector(const IndexTask &task)
      -> std::unique_ptr<Corrector> override;

 private:
  bool _forward;
  Index _start;
};

}  // namespace xfdtd

#endif  // _XFDTD_CORE_TFSF_1D_H_
