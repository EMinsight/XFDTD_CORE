#ifndef _XFDTD_LIB_TFSF_H_
#define _XFDTD_LIB_TFSF_H_

#include <xfdtd/grid_space/grid_space.h>

#include <cstddef>
#include <memory>
#include <xtensor/xtensor.hpp>

#include "xfdtd/coordinate_system/coordinate_system.h"
#include "xfdtd/waveform_source/waveform_source.h"

namespace xfdtd {

class TFSF : public WaveformSource {
 public:
  TFSF(std::size_t x, std::size_t y, std::size_t z, double theta, double phi,
       double psi, std::unique_ptr<Waveform> waveform);

  TFSF(TFSF &&) noexcept = default;

  TFSF &operator=(TFSF &&) noexcept = default;

  ~TFSF() override = default;

  void correctMaterialSpace() override;

  void correctUpdateCoefficient() override;

  void updateWaveformSourceE() override;

  void updateWaveformSourceH() override;

  std::size_t x() const;

  std::size_t y() const;

  std::size_t z() const;

  double theta() const;

  double phi() const;

  double psi() const;

  double sinTheta() const;

  double cosTheta() const;

  double sinPhi() const;

  double cosPhi() const;

  double sinPsi() const;

  double cosPsi() const;

  Vector k() const;

 protected:
  void defaultInit(std::shared_ptr<GridSpace> grid_space,
                   std::shared_ptr<CalculationParam> calculation_param,
                   std::shared_ptr<EMF> emf) override;

  const GridBox *boxPtr() const;

  double exInc(std::size_t i, std::size_t j, std::size_t k);

  double eyInc(std::size_t i, std::size_t j, std::size_t k);

  double ezInc(std::size_t i, std::size_t j, std::size_t k);

  double hxInc(std::size_t i, std::size_t j, std::size_t k);

  double hyInc(std::size_t i, std::size_t j, std::size_t k);

  double hzInc(std::size_t i, std::size_t j, std::size_t k);

  double cax();

  double cay();

  double caz();

  double cbx();

  double cby();

  double cbz();

 protected:
  xt::xarray<double> _projection_x_int;
  xt::xarray<double> _projection_y_int;
  xt::xarray<double> _projection_z_int;
  xt::xarray<double> _projection_x_half;
  xt::xarray<double> _projection_y_half;
  xt::xarray<double> _projection_z_half;

  xt::xarray<double> _ex_inc;
  xt::xarray<double> _ey_inc;
  xt::xarray<double> _ez_inc;
  xt::xarray<double> _hx_inc;
  xt::xarray<double> _hy_inc;
  xt::xarray<double> _hz_inc;

 private:
  std::size_t _x, _y, _z;
  double _theta, _phi, _psi;
  double _sin_theta, _cos_theta, _sin_phi, _cos_phi, _sin_psi, _cos_psi;
  Vector _k;
  xt::xtensor<double, 2> _rotation_matrix;
  Vector _k_e;
  xt::xarray<double> _transform_e, _transform_h;

  std::unique_ptr<GridBox> _box;
  double _ratio_delta;
  std::size_t _auxiliary_size;

  // IFA
  xt::xarray<double> _e_inc;
  xt::xarray<double> _h_inc;

  double _scaled_dl;
  double _ceie;
  double _chih;
  double _ceihi;
  double _chiei;
  double _abc_coff_0, _abc_coff_1;
  double _a, _b;

  void initTransform();

  void calculateProjection();

  Grid calculateInjectPostion();
};

}  // namespace xfdtd

#endif  // _XFDTD_LIB_TFSF_H_
