#ifndef _XFDTD_LIB_NFFFT_H_
#define _XFDTD_LIB_NFFFT_H_

#include <complex>
#include <cstddef>

#include "xfdtd/calculation_param/calculation_param.h"
#include "xfdtd/coordinate_system/coordinate_system.h"
#include "xfdtd/electromagnetic_field/electromagnetic_field.h"
#include "xfdtd/grid_space/grid_space.h"

namespace xfdtd {

class NFFFT {
 public:
  NFFFT(std::size_t distance_x, std::size_t distance_y, std::size_t distance_z,
        xt::xarray<double> frequencies, std::string output_dir);

  NFFFT(const NFFFT&) = delete;

  NFFFT(NFFFT&&) noexcept = default;

  NFFFT& operator=(const NFFFT&) = delete;

  NFFFT& operator=(NFFFT&&) noexcept = default;

  ~NFFFT() = default;

  void init(std::shared_ptr<const GridSpace> grid_space,
            std::shared_ptr<const CalculationParam> calculation_param,
            std::shared_ptr<const EMF> emf);

  void update();

  void output();

  void outputRadiationPower();

  void processFarField(const xt::xarray<double>& theta,
                       const xt::xarray<double>& phi,
                       const std::string& sub_dir,
                       const Vector& origin = Vector{0.0, 0.0, 0.0});

  void setOutputDir(const std::string& out_dir);

 protected:
  const GridSpace* gridSpacePtr() const;

  const CalculationParam* calculationParamPtr() const;

  const EMF* emfPtr() const;

  const GridBox* gridBoxPtr() const;

 private:
  std::size_t _distance_x, _distance_y, _distance_z;
  xt::xarray<double> _frequencies;
  std::string _output_dir;

  std::shared_ptr<const GridSpace> _grid_space;
  std::shared_ptr<const CalculationParam> _calculation_param;
  std::shared_ptr<const EMF> _emf;

  std::unique_ptr<GridBox> _grid_box;

  xt::xarray<std::complex<double>> _jx_yn, _jx_yp, _jx_zn, _jx_zp;
  xt::xarray<std::complex<double>> _jy_xn, _jy_xp, _jy_zn, _jy_zp;
  xt::xarray<std::complex<double>> _jz_xn, _jz_xp, _jz_yn, _jz_yp;
  xt::xarray<std::complex<double>> _mx_yn, _mx_yp, _mx_zn, _mx_zp;
  xt::xarray<std::complex<double>> _my_xn, _my_xp, _my_zn, _my_zp;
  xt::xarray<std::complex<double>> _mz_xn, _mz_xp, _mz_yn, _mz_yp;

  xt::xarray<std::complex<double>> _transform_e, _transform_h;

  xt::xarray<std::complex<double>> _a_theta, _a_phi, _f_theta, _f_phi;

  void initDFT();

  void vectorPotential(
      const xt::xarray<double>& jx_yn, const xt::xarray<double>& jx_yp,
      const xt::xarray<double>& jx_zn, const xt::xarray<double>& jx_zp,
      const xt::xarray<double>& jy_xn, const xt::xarray<double>& jy_xp,
      const xt::xarray<double>& jy_zn, const xt::xarray<double>& jy_zp,
      const xt::xarray<double>& jz_xn, const xt::xarray<double>& jz_xp,
      const xt::xarray<double>& jz_yn, const xt::xarray<double>& jz_yp,
      const xt::xarray<double>& theta, const xt::xarray<double>& phi,
      const Vector& origin);
};

}  // namespace xfdtd

#endif  // _XFDTD_LIB_NFFFT_H_
