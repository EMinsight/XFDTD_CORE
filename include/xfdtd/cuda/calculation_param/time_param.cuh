#ifndef __XFDTD_CORE_CUDA_TIME_PARAM_CUH__
#define __XFDTD_CORE_CUDA_TIME_PARAM_CUH__

#include <xfdtd/calculation_param/time_param.h>
#include <xfdtd/common/type_define.h>

#include <iostream>
#include <xfdtd/cuda/common.cuh>
#include <xfdtd/cuda/tensor.cuh>
#include <xfdtd/cuda/tensor_hd.cuh>

namespace xfdtd {

namespace cuda {

class TimeParam {
  friend class TimeParamHD;

 public:
  TimeParam() = default;

  XFDTD_CORE_CUDA_DUAL TimeParam(Real dt, Index start_time_step, Index size,
                                 Index current_time_step)
      : _dt{dt},
        _start_time_step{start_time_step},
        _size{size},
        _current_time_step{current_time_step} {}

  XFDTD_CORE_CUDA_DUAL auto dt() const { return _dt; }

  XFDTD_CORE_CUDA_DUAL auto startTimeStep() const { return _start_time_step; }

  XFDTD_CORE_CUDA_DUAL auto size() const { return _size; }

  XFDTD_CORE_CUDA_DUAL auto currentTimeStep() const {
    return _current_time_step;
  }

  XFDTD_CORE_CUDA_DUAL auto endTimeStep() const {
    return _start_time_step + _size;
  }

  XFDTD_CORE_CUDA_DUAL auto remainingTimeStep() const {
    return _start_time_step + _size - _current_time_step;
  }

  XFDTD_CORE_CUDA_DUAL auto nextTimeStep() { return ++_current_time_step; }

  /**
   * @brief [1,2,3,...,size] * dt
   */
  XFDTD_CORE_CUDA_DUAL auto eTime() const {
    auto interval = dt();
    auto e_time = Tensor<Real, 1>::from_shape({size()});
    for (Index i = 0; i < _size; ++i) {
      e_time(i) = interval * (i + 1);
    }

    return e_time;
  }

  /**
   * @brief [0.5,1.5,2.5,...,size-0.5] * dt
   */
  XFDTD_CORE_CUDA_DUAL auto hTime() const {
    auto interval = dt();
    auto h_time = Tensor<Real, 1>::from_shape({size()});
    for (Index i = 0; i < _size; ++i) {
      h_time(i) = (interval + 0.5) * i;
    }

    return h_time;
  }

 private:
  Real _dt{};
  Index _start_time_step{}, _size{}, _current_time_step{};
};

class TimeParamHD {
 public:
  using HostTimeParam = xfdtd::TimeParam;
  using DeviceTimeParam = xfdtd::cuda::TimeParam;

 public:
  TimeParamHD(HostTimeParam* time_param) : _time_param{time_param} {}

  TimeParamHD(const TimeParamHD&) = delete;

  TimeParamHD(TimeParamHD&& ohter) noexcept
      : _time_param{ohter._time_param},
        _device_time_param{ohter._device_time_param} {
    ohter._device_time_param = nullptr;
  }

  ~TimeParamHD() { releaseDevice(); }

  auto operator=(const TimeParamHD&) -> TimeParamHD& = delete;

  auto operator=(TimeParamHD&& other) noexcept -> TimeParamHD& {
    if (this != &other) {
      _time_param = std::move(other._time_param);
      _device_time_param = other._device_time_param;
      other._device_time_param = nullptr;
    }

    return *this;
  }

  auto setHostTimeParam(HostTimeParam* time_param) { _time_param = time_param; }

  auto copyHostToDevice() {
    auto d = DeviceTimeParam{};
    d._dt = _time_param->dt();
    d._start_time_step = _time_param->startTimeStep();
    d._size = _time_param->size();
    d._current_time_step = _time_param->currentTimeStep();

    if (_device_time_param != nullptr) {
      XFDTD_CORE_CUDA_CHECK_CUDA_ERROR(cudaFree(_device_time_param));
    }

    XFDTD_CORE_CUDA_CHECK_CUDA_ERROR(
        cudaMalloc(&_device_time_param, sizeof(DeviceTimeParam)));

    XFDTD_CORE_CUDA_CHECK_CUDA_ERROR(cudaMemcpy(_device_time_param, &d,
                                                sizeof(DeviceTimeParam),
                                                cudaMemcpyHostToDevice));
  }

  auto copyDeviceToHost() {
    if (_device_time_param == nullptr) {
      throw std::runtime_error("DeviceTimeParam is nullptr");
    }

    auto d = DeviceTimeParam{};

    XFDTD_CORE_CUDA_CHECK_CUDA_ERROR(cudaMemcpy(&d, _device_time_param,
                                                sizeof(DeviceTimeParam),
                                                cudaMemcpyDeviceToHost));
    // TODO: copy d to _time_param
    std::cout << "TODO: TimeParamHD::copyDeviceToHost Set Current Time Step\n";
  }

  auto moveBackToHost() -> void {}

  auto releaseDevice() -> void {
    if (_device_time_param != nullptr) {
      XFDTD_CORE_CUDA_CHECK_CUDA_ERROR(cudaFree(_device_time_param));
      _device_time_param = nullptr;
    }
  }

  auto host() { return _time_param; }

  auto device() { return _device_time_param; }

  auto host() const { return _time_param; }

  auto device() const { return _device_time_param; }


 private:
  HostTimeParam* _time_param{};
  DeviceTimeParam* _device_time_param{};
};

}  // namespace cuda

}  // namespace xfdtd

#endif  // __XFDTD_CORE_CUDA_TIME_PARAM_CUH__
