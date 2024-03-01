#ifndef _XFDTD_LIB_WAVEFORM_SOURCE_CORRECTOR_H_
#define _XFDTD_LIB_WAVEFORM_SOURCE_CORRECTOR_H_

#include <memory>
#include <utility>

#include "corrector/corrector.h"
#include "divider/divider.h"
#include "xfdtd/calculation_param/calculation_param.h"
#include "xfdtd/electromagnetic_field/electromagnetic_field.h"
#include "xfdtd/grid_space/grid_space.h"

namespace xfdtd {

class WaveformSourceCorrector : public Corrector {
 public:
  WaveformSourceCorrector(Divider::IndexTask task,
                          Divider::IndexTask local_task,
                          std::shared_ptr<GridSpace> grid_space,
                          std::shared_ptr<CalculationParam> calculation_param,
                          std::shared_ptr<EMF> emf,
                          const xt::xarray<double>& waveform)
      : _task{std::move(task)},
        _local_task{std::move(local_task)},
        _grid_space{std::move(grid_space)},
        _calculation_param{std::move(calculation_param)},
        _emf{std::move(emf)},
        _waveform{waveform} {}

  ~WaveformSourceCorrector() override = default;

  auto task() -> Divider::IndexTask { return _task; }

  auto localTask() -> Divider::IndexTask { return _local_task; }

  auto gridSpace() -> std::shared_ptr<GridSpace> { return _grid_space; }

  auto calculationParam() -> std::shared_ptr<CalculationParam> {
    return _calculation_param;
  }

  auto emf() -> std::shared_ptr<EMF> { return _emf; }

  auto waveform() -> const xt::xarray<double>& { return _waveform; }

 protected:
  auto calculationParamPtr() { return _calculation_param.get(); }

  auto emfPtr() { return _emf.get(); }

 private:
  Divider::IndexTask _task, _local_task;
  std::shared_ptr<GridSpace> _grid_space;
  std::shared_ptr<CalculationParam> _calculation_param;
  std::shared_ptr<EMF> _emf;
  const xt::xarray<double>& _waveform;
};

class TFSFCorrector : public WaveformSourceCorrector {
 public:
  TFSFCorrector(Divider::IndexTask task, Divider::IndexTask local_task,
                std::shared_ptr<GridSpace> grid_space,
                std::shared_ptr<CalculationParam> calculation_param,
                std::shared_ptr<EMF> emf, const xt::xarray<double>& waveform,
                Grid start, Divider::IndexTask xn_task,
                Divider::IndexTask xp_task, Divider::IndexTask yn_task,
                Divider::IndexTask yp_task, Divider::IndexTask zn_task,
                Divider::IndexTask zp_task,
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
      : WaveformSourceCorrector{std::move(task),
                                std::move(local_task),
                                std::move(grid_space),
                                std::move(calculation_param),
                                std::move(emf),
                                waveform},
        _start{start},
        _xn_task{std::move(xn_task)},
        _xp_task{std::move(xp_task)},
        _yn_task{std::move(yn_task)},
        _yp_task{std::move(yp_task)},
        _zn_task{std::move(zn_task)},
        _zp_task{std::move(zp_task)},
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

  Grid start() const { return _start; }

 protected:
  std::size_t jStartXN() const { return _xn_task._y_range[0]; }

  std::size_t jEndXN() const { return _xn_task._y_range[1]; }

  std::size_t kStartXN() const { return _xn_task._z_range[0]; }

  std::size_t kEndXN() const { return _xn_task._z_range[1]; }

  std::size_t jStartXP() const { return _xp_task._y_range[0]; }

  std::size_t jEndXP() const { return _xp_task._y_range[1]; }

  std::size_t kStartXP() const { return _xp_task._z_range[0]; }

  std::size_t kEndXP() const { return _xp_task._z_range[1]; }

  std::size_t iStartYN() const { return _yn_task._x_range[0]; }

  std::size_t iEndYN() const { return _yn_task._x_range[1]; }

  std::size_t kStartYN() const { return _yn_task._z_range[0]; }

  std::size_t kEndYN() const { return _yn_task._z_range[1]; }

  std::size_t iStartYP() const { return _yp_task._x_range[0]; }

  std::size_t iEndYP() const { return _yp_task._x_range[1]; }

  std::size_t kStartYP() const { return _yp_task._z_range[0]; }

  std::size_t kEndYP() const { return _yp_task._z_range[1]; }

  std::size_t iStartZN() const { return _zn_task._x_range[0]; }

  std::size_t iEndZN() const { return _zn_task._x_range[1]; }

  std::size_t jStartZN() const { return _zn_task._y_range[0]; }

  std::size_t jEndZN() const { return _zn_task._y_range[1]; }

  std::size_t iStartZP() const { return _zp_task._x_range[0]; }

  std::size_t iEndZP() const { return _zp_task._x_range[1]; }

  std::size_t jStartZP() const { return _zp_task._y_range[0]; }

  std::size_t jEndZP() const { return _zp_task._y_range[1]; }

 private:
  Grid _start;
  Divider::IndexTask _xn_task, _xp_task, _yn_task, _yp_task, _zn_task, _zp_task;

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
  TFSF1DCorrector(Divider::IndexTask task, Divider::IndexTask local_task,
                  std::shared_ptr<GridSpace> grid_space,
                  std::shared_ptr<CalculationParam> calculation_param,
                  std::shared_ptr<EMF> emf, const xt::xarray<double>& waveform,
                  Grid start, Divider::IndexTask xn_task,
                  Divider::IndexTask xp_task, Divider::IndexTask yn_task,
                  Divider::IndexTask yp_task, Divider::IndexTask zn_task,
                  Divider::IndexTask zp_task,
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
      : TFSFCorrector{std::move(task),
                      std::move(local_task),
                      std::move(grid_space),
                      std::move(calculation_param),
                      std::move(emf),
                      waveform,
                      start,
                      std::move(xn_task),
                      std::move(xp_task),
                      std::move(yn_task),
                      std::move(yp_task),
                      std::move(zn_task),
                      std::move(zp_task),
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

  ~TFSF1DCorrector() override = default;

  void correctE() override;

  void correctH() override;

 private:
};

class TFSF2DCorrector : public TFSFCorrector {
 public:
  TFSF2DCorrector(Divider::IndexTask task, Divider::IndexTask local_task,
                  std::shared_ptr<GridSpace> grid_space,
                  std::shared_ptr<CalculationParam> calculation_param,
                  std::shared_ptr<EMF> emf, const xt::xarray<double>& waveform,
                  Grid start, Divider::IndexTask xn_task,
                  Divider::IndexTask xp_task, Divider::IndexTask yn_task,
                  Divider::IndexTask yp_task, Divider::IndexTask zn_task,
                  Divider::IndexTask zp_task,
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
      : TFSFCorrector{std::move(task),
                      std::move(local_task),
                      std::move(grid_space),
                      std::move(calculation_param),
                      std::move(emf),
                      waveform,
                      start,
                      std::move(xn_task),
                      std::move(xp_task),
                      std::move(yn_task),
                      std::move(yp_task),
                      std::move(zn_task),
                      std::move(zp_task),
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
  TFSF3DCorrector(Divider::IndexTask task, Divider::IndexTask local_task,
                  std::shared_ptr<GridSpace> grid_space,
                  std::shared_ptr<CalculationParam> calculation_param,
                  std::shared_ptr<EMF> emf, const xt::xarray<double>& waveform,
                  Grid start, Divider::IndexTask xn_task,
                  Divider::IndexTask xp_task, Divider::IndexTask yn_task,
                  Divider::IndexTask yp_task, Divider::IndexTask zn_task,
                  Divider::IndexTask zp_task,
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
      : TFSFCorrector{std::move(task),
                      std::move(local_task),
                      std::move(grid_space),
                      std::move(calculation_param),
                      std::move(emf),
                      waveform,
                      start,
                      std::move(xn_task),
                      std::move(xp_task),
                      std::move(yn_task),
                      std::move(yp_task),
                      std::move(zn_task),
                      std::move(zp_task),
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

#endif  // _XFDTD_LIB_WAVEFORM_SOURCE_CORRECTOR_H_
