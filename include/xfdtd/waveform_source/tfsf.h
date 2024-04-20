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

  void updateWaveformSource() override;

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

  IndexTask taskXN() const;

  IndexTask taskXP() const;

  IndexTask taskYN() const;

  IndexTask taskYP() const;

  IndexTask taskZN() const;

  IndexTask taskZP() const;

  IndexTask globalEyTaskXN() const;

  IndexTask globalEzTaskXN() const;

  IndexTask globalEyTaskXP() const;

  IndexTask globalEzTaskXP() const;

  IndexTask globalEzTaskYN() const;

  IndexTask globalExTaskYN() const;

  IndexTask globalEzTaskYP() const;

  IndexTask globalExTaskYP() const;

  IndexTask globalExTaskZN() const;

  IndexTask globalEyTaskZN() const;

  IndexTask globalExTaskZP() const;

  IndexTask globalEyTaskZP() const;

 protected:
  void defaultInit(std::shared_ptr<GridSpace> grid_space,
                   std::shared_ptr<CalculationParam> calculation_param,
                   std::shared_ptr<EMF> emf) override;

  Real exInc(std::size_t i, std::size_t j, std::size_t k);

  Real eyInc(std::size_t i, std::size_t j, std::size_t k);

  Real ezInc(std::size_t i, std::size_t j, std::size_t k);

  Real hxInc(std::size_t i, std::size_t j, std::size_t k);

  Real hyInc(std::size_t i, std::size_t j, std::size_t k);

  Real hzInc(std::size_t i, std::size_t j, std::size_t k);

  Real cax();

  Real cay();

  Real caz();

  Real cbx();

  Real cby();

  Real cbz();

  IndexTask nodeEyTaskXN(const IndexTask &task) const;

  IndexTask nodeEzTaskXN(const IndexTask &task) const;

  IndexTask nodeEyTaskXP(const IndexTask &task) const;

  IndexTask nodeEzTaskXP(const IndexTask &task) const;

  IndexTask nodeExTaskYN(const IndexTask &task) const;

  IndexTask nodeEzTaskYN(const IndexTask &task) const;

  IndexTask nodeExTaskYP(const IndexTask &task) const;

  IndexTask nodeEzTaskYP(const IndexTask &task) const;

  IndexTask nodeExTaskZN(const IndexTask &task) const;

  IndexTask nodeEyTaskZN(const IndexTask &task) const;

  IndexTask nodeExTaskZP(const IndexTask &task) const;

  IndexTask nodeEyTaskZP(const IndexTask &task) const;

 protected:
  Array1D<Real> _projection_x_int;
  Array1D<Real> _projection_y_int;
  Array1D<Real> _projection_z_int;
  Array1D<Real> _projection_x_half;
  Array1D<Real> _projection_y_half;
  Array1D<Real> _projection_z_half;

  Array1D<Real> _ex_inc;
  Array1D<Real> _ey_inc;
  Array1D<Real> _ez_inc;
  Array1D<Real> _hx_inc;
  Array1D<Real> _hy_inc;
  Array1D<Real> _hz_inc;

 private:
  std::size_t _x, _y, _z;
  Real _theta, _phi, _psi;
  Real _sin_theta, _cos_theta, _sin_phi, _cos_phi, _sin_psi, _cos_psi;
  Vector _k;
  Array2D<Real> _rotation_matrix;
  Vector _k_e;
  Array<Real> _transform_e, _transform_h;

  //   std::unique_ptr<GridBox> _box;
  GridBox _global_box;
  Real _ratio_delta;
  std::size_t _auxiliary_size;

  // IFA
  Array1D<Real> _e_inc;
  Array1D<Real> _h_inc;

  Real _scaled_dl;
  Real _ceie;
  Real _chih;
  Real _ceihi;
  Real _chiei;
  Real _abc_coff_0, _abc_coff_1;
  Real _a{}, _b{};  // !!!!!!!!! Don't forget to initialize these variables

  void initTransform();

  void calculateProjection();

  Grid calculateInjectPostion();
};

}  // namespace xfdtd

#endif  // _XFDTD_CORE_TFSF_H_
