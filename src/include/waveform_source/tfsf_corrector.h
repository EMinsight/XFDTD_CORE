#ifndef __XFDTD_CORE_TFSF_CORRECTOR_H__
#define __XFDTD_CORE_TFSF_CORRECTOR_H__

#include <xfdtd/common/index_task.h>
#include <xfdtd/common/type_define.h>

#include "waveform_source/waveform_source_corrector.h"

namespace xfdtd {

class TFSFCorrector : public WaveformSourceCorrector {
 public:
  TFSFCorrector(Grid global_start, std::shared_ptr<GridSpace> grid_space,
                std::shared_ptr<CalculationParam> calculation_param,
                std::shared_ptr<EMF> emf, const Array1D<Real>& waveform,
                IndexTask global_ey_task_xn, IndexTask global_ez_task_xn,
                IndexTask global_ey_task_xp, IndexTask global_ez_task_xp,
                IndexTask global_ez_task_yn, IndexTask global_ex_task_yn,
                IndexTask global_ez_task_yp, IndexTask global_ex_task_yp,
                IndexTask global_ex_task_zn, IndexTask global_ey_task_zn,
                IndexTask global_ex_task_zp, IndexTask global_ey_task_zp,
                const Array1D<Real>& projection_x_int,
                const Array1D<Real>& projection_y_int,
                const Array1D<Real>& projection_z_int,
                const Array1D<Real>& projection_x_half,
                const Array1D<Real>& projection_y_half,
                const Array1D<Real>& projection_z_half, Array1D<Real>& ex_inc,
                Array1D<Real>& ey_inc, Array1D<Real>& ez_inc,
                Array1D<Real>& hx_inc, Array1D<Real>& hy_inc,
                Array1D<Real>& hz_inc, Real cax, Real cbx, Real cay, Real cby,
                Real caz, Real cbz)
      : WaveformSourceCorrector{IndexTask{},
                                IndexTask{},
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

  Real exInc(std::size_t i, std::size_t j, std::size_t k) const;

  Real eyInc(std::size_t i, std::size_t j, std::size_t k) const;

  Real ezInc(std::size_t i, std::size_t j, std::size_t k) const;

  Real hxInc(std::size_t i, std::size_t j, std::size_t k) const;

  Real hyInc(std::size_t i, std::size_t j, std::size_t k) const;

  Real hzInc(std::size_t i, std::size_t j, std::size_t k) const;

  Real cax() const { return _cax; }

  Real cay() const { return _cay; }

  Real caz() const { return _caz; }

  Real cbx() const { return _cbx; }

  Real cby() const { return _cby; }

  Real cbz() const { return _cbz; }

  Grid globalStart() const { return _global_start; }

  std::string toString() const override;

  IndexTask globalEyTaskXN() const;

  IndexTask globalEzTaskXN() const;

  IndexTask globalEyTaskXP() const;

  IndexTask globalEzTaskXP() const;

  IndexTask globalExTaskYN() const;

  IndexTask globalEzTaskYN() const;

  IndexTask globalExTaskYP() const;

  IndexTask globalEzTaskYP() const;

  IndexTask globalExTaskZN() const;

  IndexTask globalEyTaskZN() const;

  IndexTask globalExTaskZP() const;

  IndexTask globalEyTaskZP() const;

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
  IndexTask _global_ey_task_xn, _global_ez_task_xn, _global_ey_task_xp,
      _global_ez_task_xp, _global_ez_task_yn, _global_ex_task_yn,
      _global_ez_task_yp, _global_ex_task_yp, _global_ex_task_zn,
      _global_ey_task_zn, _global_ex_task_zp, _global_ey_task_zp;

  const Array1D<Real>& _projection_x_int;
  const Array1D<Real>& _projection_y_int;
  const Array1D<Real>& _projection_z_int;
  const Array1D<Real>& _projection_x_half;
  const Array1D<Real>& _projection_y_half;
  const Array1D<Real>& _projection_z_half;

  // IFA
  Array1D<Real>& _ex_inc;
  Array1D<Real>& _ey_inc;
  Array1D<Real>& _ez_inc;
  Array1D<Real>& _hx_inc;
  Array1D<Real>& _hy_inc;
  Array1D<Real>& _hz_inc;
  const Real _cax, _cbx, _cay, _cby, _caz, _cbz;

  Real _a, _b;
};

class TFSF1DCorrector : public TFSFCorrector {
 public:
  TFSF1DCorrector(Grid global_start, std::shared_ptr<GridSpace> grid_space,
                  std::shared_ptr<CalculationParam> calculation_param,
                  std::shared_ptr<EMF> emf, const Array1D<Real>& waveform,
                  bool forward, const Array1D<Real>& projection_x_int,
                  const Array1D<Real>& projection_y_int,
                  const Array1D<Real>& projection_z_int,
                  const Array1D<Real>& projection_x_half,
                  const Array1D<Real>& projection_y_half,
                  const Array1D<Real>& projection_z_half, Array1D<Real>& ex_inc,
                  Array1D<Real>& ey_inc, Array1D<Real>& ez_inc,
                  Array1D<Real>& hx_inc, Array1D<Real>& hy_inc,
                  Array1D<Real>& hz_inc, Real cax, Real cbx, Real cay, Real cby,
                  Real caz, Real cbz)
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
                      makeIndexTask(makeIndexRange(0, 1), makeIndexRange(0, 1),
                                    makeIndexRange(global_start.k(),
                                                   global_start.k() + 1)),
                      {},
                      makeIndexTask(makeIndexRange(0, 1), makeIndexRange(0, 1),
                                    makeIndexRange(global_start.k(),
                                                   global_start.k() + 1)),
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
                      cbz},
        _forward{forward} {};

  ~TFSF1DCorrector() override = default;

  void correctE() override;

  void correctH() override;

 private:
  const bool _forward;
};

class TFSF2DCorrector : public TFSFCorrector {
 public:
  TFSF2DCorrector(Grid global_start, std::shared_ptr<GridSpace> grid_space,
                  std::shared_ptr<CalculationParam> calculation_param,
                  std::shared_ptr<EMF> emf, const Array1D<Real>& waveform,
                  IndexTask global_ez_task_xn, IndexTask global_ez_task_xp,
                  IndexTask global_ez_task_yn, IndexTask global_ez_task_yp,
                  const Array1D<Real>& projection_x_int,
                  const Array1D<Real>& projection_y_int,
                  const Array1D<Real>& projection_z_int,
                  const Array1D<Real>& projection_x_half,
                  const Array1D<Real>& projection_y_half,
                  const Array1D<Real>& projection_z_half, Array1D<Real>& ex_inc,
                  Array1D<Real>& ey_inc, Array1D<Real>& ez_inc,
                  Array1D<Real>& hx_inc, Array1D<Real>& hy_inc,
                  Array1D<Real>& hz_inc, Real cax, Real cbx, Real cay, Real cby,
                  Real caz, Real cbz)
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
                  std::shared_ptr<EMF> emf, const Array1D<Real>& waveform,
                  IndexTask global_ey_task_xn, IndexTask global_ez_task_xn,
                  IndexTask global_ey_task_xp, IndexTask global_ez_task_xp,
                  IndexTask global_ez_task_yn, IndexTask global_ex_task_yn,
                  IndexTask global_ez_task_yp, IndexTask global_ex_task_yp,
                  IndexTask global_ex_task_zn, IndexTask global_ey_task_zn,
                  IndexTask global_ex_task_zp, IndexTask global_ey_task_zp,
                  const Array1D<Real>& projection_x_int,
                  const Array1D<Real>& projection_y_int,
                  const Array1D<Real>& projection_z_int,
                  const Array1D<Real>& projection_x_half,
                  const Array1D<Real>& projection_y_half,
                  const Array1D<Real>& projection_z_half, Array1D<Real>& ex_inc,
                  Array1D<Real>& ey_inc, Array1D<Real>& ez_inc,
                  Array1D<Real>& hx_inc, Array1D<Real>& hy_inc,
                  Array1D<Real>& hz_inc, Real cax, Real cbx, Real cay, Real cby,
                  Real caz, Real cbz)
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
