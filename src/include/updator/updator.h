#ifndef _XFDTD_LIB_UPDATOR_H_
#define _XFDTD_LIB_UPDATOR_H_

#include <memory>

#include "divider/divider.h"
#include "xfdtd/calculation_param/calculation_param.h"
#include "xfdtd/electromagnetic_field/electromagnetic_field.h"
#include "xfdtd/grid_space/grid_space.h"

namespace xfdtd {

class Updator {
 public:
  Updator(std::shared_ptr<const GridSpace> grid_space,
          std::shared_ptr<const CalculationParam> calculation_param,
          std::shared_ptr<EMF> emf, Divider::IndexTask task);

  Updator(const Updator&) = default;

  Updator(Updator&&) noexcept = default;

  Updator& operator=(const Updator&) = default;

  Updator& operator=(Updator&&) noexcept = default;

  virtual ~Updator() = default;

  virtual void updateE() = 0;

  virtual void updateH() = 0;

  auto task() const { return _task; }

  virtual std::string toString() const;

  // range should be {is:ie, js-1, ks:ke}
  void setHxBufferForEdgeYN(xt::xarray<double> hx_buffer_yn) {
    _hx_buffer_yn = std::move(hx_buffer_yn);
  }

  // shape should be {is:ie, js:je, 1}
  void setHxBufferForEdgeZN(xt::xarray<double> hx_buffer_zn) {
    _hx_buffer_zn = std::move(hx_buffer_zn);
  }

  // shape should be {1, js:je, ks:ke}
  void setHyBufferForEdgeXN(xt::xarray<double> hy_buffer_xn) {
    _hy_buffer_xn = std::move(hy_buffer_xn);
  }

  // range should be {is:ie, js:je, ks-1}
  void setHyBufferForEdgeZN(xt::xarray<double> hy_buffer_zn) {
    _hy_buffer_zn = std::move(hy_buffer_zn);
  }

  // shape should be {1, js:je, ks:ke}
  void setHzBufferForEdgeXN(xt::xarray<double> hz_buffer_xn) {
    _hz_buffer_xn = std::move(hz_buffer_xn);
  }

  // shape should be {is:ie, 1, ks:ke}
  void setHzBufferForEdgeYN(xt::xarray<double> hz_buffer_yn) {
    _hz_buffer_yn = std::move(hz_buffer_yn);
  }

  // only work for shared memory model
  bool containXNEdge() const { return task().xRange().start() == 0; }

  // only work for shared memory model
  bool containYNEdge() const { return task().yRange().start() == 0; }

  // only work for shared memory model
  bool containZNEdge() const { return task().zRange().start() == 0; }

 protected:
  std::shared_ptr<const GridSpace> _grid_space;
  std::shared_ptr<const CalculationParam> _calculation_param;
  std::shared_ptr<EMF> _emf;

  virtual void updateEEdge() = 0;

  virtual void updateHEdge() = 0;

  auto&& hxBufferYN() { return _hx_buffer_yn; }

  auto&& hxBufferZN() { return _hx_buffer_zn; }

  auto&& hyBufferXN() { return _hy_buffer_xn; }

  auto&& hyBufferZN() { return _hy_buffer_zn; }

  auto&& hzBufferXN() { return _hz_buffer_xn; }

  auto&& hzBufferYN() { return _hz_buffer_yn; }

 private:
  Divider::IndexTask _task;

  xt::xarray<double> _hx_buffer_yn, _hx_buffer_zn, _hy_buffer_xn, _hy_buffer_zn,
      _hz_buffer_xn, _hz_buffer_yn;
};

}  // namespace xfdtd

#endif  // _XFDTD_LIB_UPDATOR_H_
