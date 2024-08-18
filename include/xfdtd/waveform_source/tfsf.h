#ifndef _XFDTD_CORE_TFSF_H_
#define _XFDTD_CORE_TFSF_H_

#include <xfdtd/common/index_task.h>
#include <xfdtd/common/type_define.h>
#include <xfdtd/coordinate_system/coordinate_system.h>
#include <xfdtd/grid_space/grid_space.h>
#include <xfdtd/waveform_source/waveform_source.h>

#include <cstddef>
#include <memory>

namespace xfdtd {

class TFSF : public WaveformSource {
 public:
  TFSF(std::size_t x, std::size_t y, std::size_t z, Real theta, Real phi,
       Real psi, std::unique_ptr<Waveform> waveform);

  TFSF(TFSF &&) noexcept = default;

  TFSF &operator=(TFSF &&) noexcept = default;

  ~TFSF() override = default;

  void correctMaterialSpace() override;

  void correctUpdateCoefficient() override;

  void initTimeDependentVariable() override;

  std::size_t x() const;

  std::size_t y() const;

  std::size_t z() const;

  Real theta() const;

  Real phi() const;

  Real psi() const;

  Real sinTheta() const;

  Real cosTheta() const;

  Real sinPhi() const;

  Real cosPhi() const;

  Real sinPsi() const;

  Real cosPsi() const;

  Vector k() const;

  GridBox globalBox() const;

  auto globalTask() const -> IndexTask;

  Real cax() const;

  Real cay() const;

  Real caz() const;

  Real cbx() const;

  Real cby() const;

  Real cbz() const;

  auto transformE() const -> Vector { return _transform_e; }

  auto transformH() const -> Vector { return _transform_h; }

  auto projectionXInt() const -> const Array1D<Real> & {
    return _projection_x_int;
  }

  auto projectionYInt() const -> const Array1D<Real> & {
    return _projection_y_int;
  }

  auto projectionZInt() const -> const Array1D<Real> & {
    return _projection_z_int;
  }

  auto projectionXHalf() const -> const Array1D<Real> & {
    return _projection_x_half;
  }

  auto projectionYHalf() const -> const Array1D<Real> & {
    return _projection_y_half;
  }

  auto projectionZHalf() const -> const Array1D<Real> & {
    return _projection_z_half;
  }

  auto eInc() const -> const Array2D<Real> & { return _e_inc; }

  auto hInc() const -> const Array2D<Real> & { return _h_inc; }

  auto nodeTask() const -> IndexTask;

 protected:
  void defaultInit(std::shared_ptr<GridSpace> grid_space,
                   std::shared_ptr<CalculationParam> calculation_param,
                   std::shared_ptr<EMF> emf) override;

  auto nodeGlobalTask() const -> IndexTask;

 protected:
  Array1D<Real> _projection_x_int;
  Array1D<Real> _projection_y_int;
  Array1D<Real> _projection_z_int;
  Array1D<Real> _projection_x_half;
  Array1D<Real> _projection_y_half;
  Array1D<Real> _projection_z_half;

  // IFA
  Array2D<Real> _e_inc, _h_inc;

  Vector _transform_e, _transform_h;

 private:
  Index _x, _y, _z;
  Real _theta, _phi, _psi;
  Real _sin_theta, _cos_theta, _sin_phi, _cos_phi, _sin_psi, _cos_psi;
  Vector _k;
  Vector _k_e;

  GridBox _global_box;
  Real _ratio_delta;
  Index _auxiliary_size;

  Real _scaled_dl;
  Real _ceie;
  Real _chih;
  Real _ceihi;
  Real _chiei;
  Real _abc_coff_0, _abc_coff_1;
  Real _a = 0, _b = 0;

  void initTransform();

  void calculateProjection();

  Grid calculateInjectPostion();
};

}  // namespace xfdtd

#endif  // _XFDTD_CORE_TFSF_H_
