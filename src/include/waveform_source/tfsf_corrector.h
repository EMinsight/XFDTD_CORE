#ifndef __XFDTD_CORE_TFSF_CORRECTOR_H__
#define __XFDTD_CORE_TFSF_CORRECTOR_H__

#include "waveform_source/waveform_source_corrector.h"
#include "xfdtd/divider/divider.h"

namespace xfdtd {

class TFSFCorrector : public WaveformSourceCorrector {
 public:
  TFSFCorrector(Grid global_start, std::shared_ptr<GridSpace> grid_space,
                std::shared_ptr<CalculationParam> calculation_param,
                std::shared_ptr<EMF> emf, const xt::xarray<double>& waveform,
                Divider::IndexTask global_ey_task_xn,
                Divider::IndexTask global_ez_task_xn,
                Divider::IndexTask global_ey_task_xp,
                Divider::IndexTask global_ez_task_xp,
                Divider::IndexTask global_ez_task_yn,
                Divider::IndexTask global_ex_task_yn,
                Divider::IndexTask global_ez_task_yp,
                Divider::IndexTask global_ex_task_yp,
                Divider::IndexTask global_ex_task_zn,
                Divider::IndexTask global_ey_task_zn,
                Divider::IndexTask global_ex_task_zp,
                Divider::IndexTask global_ey_task_zp,
                const xt::xarray<double>& projection_x_int,
                const xt::xarray<double>& projection_y_int,
                const xt::xarray<double>& projection_z_int,
                const xt::xarray<double>& projection_x_half,
                const xt::xarray<double>& projection_y_half,
                const xt::xarray<double>& projection_z_half,
                xt::xarray<double>& ex_inc, xt::xarray<double>& ey_inc,
                xt::xarray<double>& ez_inc, xt::xarray<double>& hx_inc,
                xt::xarray<double>& hy_inc, xt::xarray<double>& hz_inc,
                double cax, double cbx, double cay, double cby, double caz,
                double cbz)
      : WaveformSourceCorrector{Divider::IndexTask{},
                                Divider::IndexTask{},
                                std::move(grid_space),
                                std::move(calculation_param),
                                std::move(emf),
                                waveform},
        _global_start{global_start},
        _global_ey_task_xn{global_ey_task_xn},
        _global_ez_task_xn{global_ez_task_xn},
        _global_ey_task_xp{global_ey_task_xp},
        _global_ez_task_xp{global_ez_task_xp},
        _global_ez_task_yn{global_ez_task_yn},
        _global_ex_task_yn{global_ex_task_yn},
        _global_ez_task_yp{global_ez_task_yp},
        _global_ex_task_yp{global_ex_task_yp},
        _global_ex_task_zn{global_ex_task_zn},
        _global_ey_task_zn{global_ey_task_zn},
        _global_ex_task_zp{global_ex_task_zp},
        _global_ey_task_zp{global_ey_task_zp},
        _projection_x_int{projection_x_int},
        _projection_y_int{projection_y_int},
        _projection_z_int{projection_z_int},
        _projection_x_half{projection_x_half},
        _projection_y_half{projection_y_half},
        _projection_z_half{projection_z_half},
        _ex_inc{ex_inc},
        _ey_inc{ey_inc},
        _ez_inc{ez_inc},
        _hx_inc{hx_inc},
        _hy_inc{hy_inc},
        _hz_inc{hz_inc},
        _cax{cax},
        _cbx{cbx},
        _cay{cay},
        _cby{cby},
        _caz{caz},
        _cbz{cbz} {}

  ~TFSFCorrector() override = default;

  double exInc(std::size_t i, std::size_t j, std::size_t k) const;

  double eyInc(std::size_t i, std::size_t j, std::size_t k) const;

  double ezInc(std::size_t i, std::size_t j, std::size_t k) const;

  double hxInc(std::size_t i, std::size_t j, std::size_t k) const;

  double hyInc(std::size_t i, std::size_t j, std::size_t k) const;

  double hzInc(std::size_t i, std::size_t j, std::size_t k) const;

  double cax() const { return _cax; }

  double cay() const { return _cay; }

  double caz() const { return _caz; }

  double cbx() const { return _cbx; }

  double cby() const { return _cby; }

  double cbz() const { return _cbz; }

  Grid globalStart() const { return _global_start; }

  std::string toString() const override;

  Divider::IndexTask globalEyTaskXN() const;

  Divider::IndexTask globalEzTaskXN() const;

  Divider::IndexTask globalEyTaskXP() const;

  Divider::IndexTask globalEzTaskXP() const;

  Divider::IndexTask globalExTaskYN() const;

  Divider::IndexTask globalEzTaskYN() const;

  Divider::IndexTask globalExTaskYP() const;

  Divider::IndexTask globalEzTaskYP() const;

  Divider::IndexTask globalExTaskZN() const;

  Divider::IndexTask globalEyTaskZN() const;

  Divider::IndexTask globalExTaskZP() const;

  Divider::IndexTask globalEyTaskZP() const;

 protected:
  void correctEyXN();

  void correctEzXN();

  void correctEyXP();

  void correctEzXP();

  void correctExYN();

  void correctEzYN();

  void correctExYP();

  void correctEzYP();

  void correctExZN();

  void correctEyZN();

  void correctExZP();

  void correctEyZP();

  void correctHzXN();

  void correctHyXN();

  void correctHzXP();

  void correctHyXP();

  void correctHxYN();

  void correctHzYN();

  void correctHxYP();

  void correctHzYP();

  void correctHxZN();

  void correctHyZN();

  void correctHxZP();

  void correctHyZP();

 private:
  Grid _global_start;
  Divider::IndexTask _global_ey_task_xn, _global_ez_task_xn, _global_ey_task_xp,
      _global_ez_task_xp, _global_ez_task_yn, _global_ex_task_yn,
      _global_ez_task_yp, _global_ex_task_yp, _global_ex_task_zn,
      _global_ey_task_zn, _global_ex_task_zp, _global_ey_task_zp;

  const xt::xarray<double>& _projection_x_int;
  const xt::xarray<double>& _projection_y_int;
  const xt::xarray<double>& _projection_z_int;
  const xt::xarray<double>& _projection_x_half;
  const xt::xarray<double>& _projection_y_half;
  const xt::xarray<double>& _projection_z_half;

  // IFA
  xt::xarray<double>& _ex_inc;
  xt::xarray<double>& _ey_inc;
  xt::xarray<double>& _ez_inc;
  xt::xarray<double>& _hx_inc;
  xt::xarray<double>& _hy_inc;
  xt::xarray<double>& _hz_inc;
  const double _cax, _cbx, _cay, _cby, _caz, _cbz;

  double _a, _b;
};

class TFSF1DCorrector : public TFSFCorrector {
 public:
  TFSF1DCorrector(Grid global_start, std::shared_ptr<GridSpace> grid_space,
                  std::shared_ptr<CalculationParam> calculation_param,
                  std::shared_ptr<EMF> emf, const xt::xarray<double>& waveform,
                  Divider::IndexTask zn_task, Divider::IndexTask zp_task,
                  const xt::xarray<double>& projection_x_int,
                  const xt::xarray<double>& projection_y_int,
                  const xt::xarray<double>& projection_z_int,
                  const xt::xarray<double>& projection_x_half,
                  const xt::xarray<double>& projection_y_half,
                  const xt::xarray<double>& projection_z_half,
                  xt::xarray<double>& ex_inc, xt::xarray<double>& ey_inc,
                  xt::xarray<double>& ez_inc, xt::xarray<double>& hx_inc,
                  xt::xarray<double>& hy_inc, xt::xarray<double>& hz_inc,
                  double cax, double cbx, double cay, double cby, double caz,
                  double cbz)
      : TFSFCorrector{global_start,
                      std::move(grid_space),
                      std::move(calculation_param),
                      std::move(emf),
                      waveform,
                      {},
                      {},
                      {},
                      {},
                      {},
                      {},
                      {},
                      {},
                      {},
                      {},
                      {},
                      {},
                      projection_x_int,
                      projection_y_int,
                      projection_z_int,
                      projection_x_half,
                      projection_y_half,
                      projection_z_half,
                      ex_inc,
                      ey_inc,
                      ez_inc,
                      hx_inc,
                      hy_inc,
                      hz_inc,
                      cax,
                      cbx,
                      cay,
                      cby,
                      caz,
                      cbz} {
    throw std::runtime_error("Not implemented");
  };

  ~TFSF1DCorrector() override = default;

  void correctE() override;

  void correctH() override;

 private:
};

class TFSF2DCorrector : public TFSFCorrector {
 public:
  TFSF2DCorrector(Grid global_start, std::shared_ptr<GridSpace> grid_space,
                  std::shared_ptr<CalculationParam> calculation_param,
                  std::shared_ptr<EMF> emf, const xt::xarray<double>& waveform,
                  Divider::IndexTask global_ez_task_xn,
                  Divider::IndexTask global_ez_task_xp,
                  Divider::IndexTask global_ez_task_yn,
                  Divider::IndexTask global_ez_task_yp,
                  const xt::xarray<double>& projection_x_int,
                  const xt::xarray<double>& projection_y_int,
                  const xt::xarray<double>& projection_z_int,
                  const xt::xarray<double>& projection_x_half,
                  const xt::xarray<double>& projection_y_half,
                  const xt::xarray<double>& projection_z_half,
                  xt::xarray<double>& ex_inc, xt::xarray<double>& ey_inc,
                  xt::xarray<double>& ez_inc, xt::xarray<double>& hx_inc,
                  xt::xarray<double>& hy_inc, xt::xarray<double>& hz_inc,
                  double cax, double cbx, double cay, double cby, double caz,
                  double cbz)
      : TFSFCorrector{global_start,
                      std::move(grid_space),
                      std::move(calculation_param),
                      std::move(emf),
                      waveform,
                      {},
                      global_ez_task_xn,
                      {},
                      global_ez_task_xp,
                      global_ez_task_yn,
                      {},
                      global_ez_task_yp,
                      {},
                      {},
                      {},
                      {},
                      {},
                      projection_x_int,
                      projection_y_int,
                      projection_z_int,
                      projection_x_half,
                      projection_y_half,
                      projection_z_half,
                      ex_inc,
                      ey_inc,
                      ez_inc,
                      hx_inc,
                      hy_inc,
                      hz_inc,
                      cax,
                      cbx,
                      cay,
                      cby,
                      caz,
                      cbz} {};

  ~TFSF2DCorrector() override = default;

  void correctE() override;

  void correctH() override;

 private:
};

class TFSF3DCorrector : public TFSFCorrector {
 public:
  TFSF3DCorrector(Grid global_start, std::shared_ptr<GridSpace> grid_space,
                  std::shared_ptr<CalculationParam> calculation_param,
                  std::shared_ptr<EMF> emf, const xt::xarray<double>& waveform,
                  Divider::IndexTask global_ey_task_xn,
                  Divider::IndexTask global_ez_task_xn,
                  Divider::IndexTask global_ey_task_xp,
                  Divider::IndexTask global_ez_task_xp,
                  Divider::IndexTask global_ez_task_yn,
                  Divider::IndexTask global_ex_task_yn,
                  Divider::IndexTask global_ez_task_yp,
                  Divider::IndexTask global_ex_task_yp,
                  Divider::IndexTask global_ex_task_zn,
                  Divider::IndexTask global_ey_task_zn,
                  Divider::IndexTask global_ex_task_zp,
                  Divider::IndexTask global_ey_task_zp,
                  const xt::xarray<double>& projection_x_int,
                  const xt::xarray<double>& projection_y_int,
                  const xt::xarray<double>& projection_z_int,
                  const xt::xarray<double>& projection_x_half,
                  const xt::xarray<double>& projection_y_half,
                  const xt::xarray<double>& projection_z_half,
                  xt::xarray<double>& ex_inc, xt::xarray<double>& ey_inc,
                  xt::xarray<double>& ez_inc, xt::xarray<double>& hx_inc,
                  xt::xarray<double>& hy_inc, xt::xarray<double>& hz_inc,
                  double cax, double cbx, double cay, double cby, double caz,
                  double cbz)
      : TFSFCorrector{global_start,
                      std::move(grid_space),
                      std::move(calculation_param),
                      std::move(emf),
                      waveform,
                      global_ey_task_xn,
                      global_ez_task_xn,
                      global_ey_task_xp,
                      global_ez_task_xp,
                      global_ez_task_yn,
                      global_ex_task_yn,
                      global_ez_task_yp,
                      global_ex_task_yp,
                      global_ex_task_zn,
                      global_ey_task_zn,
                      global_ex_task_zp,
                      global_ey_task_zp,
                      projection_x_int,
                      projection_y_int,
                      projection_z_int,
                      projection_x_half,
                      projection_y_half,
                      projection_z_half,
                      ex_inc,
                      ey_inc,
                      ez_inc,
                      hx_inc,
                      hy_inc,
                      hz_inc,
                      cax,
                      cbx,
                      cay,
                      cby,
                      caz,
                      cbz} {};

  ~TFSF3DCorrector() override = default;

  void correctE() override;

  void correctH() override;

 private:
};

}  // namespace xfdtd

#endif  // __XFDTD_CORE_TFSF_CORRECTOR_H__
