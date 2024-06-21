#ifndef __XFDTD_CORE_CUDA_TENSOR_HD_CUH__
#define __XFDTD_CORE_CUDA_TENSOR_HD_CUH__

#include <xfdtd/exception/exception.h>

#include <cassert>
#include <iostream>
#include <sstream>
#include <string>
#include <xfdtd/cuda/common.cuh>
#include <xfdtd/cuda/fixed_array.cuh>
#include <xfdtd/cuda/memory.cuh>
#include <xfdtd/cuda/tensor.cuh>

namespace xfdtd {

namespace cuda {

// define marco for check cuda error
#define XFDTD_CORE_CUDA_CHECK_CUDA_ERROR(call)                             \
  {                                                                        \
    auto err = call;                                                       \
    if (err != cudaSuccess) {                                              \
      std::stringstream ss;                                                \
      ss << "CUDA error in file '" << __FILE__ << "' in line " << __LINE__ \
         << ": " << cudaGetErrorString(err) << "\n";                       \
      std::cerr << ss.str();                                               \
      throw std::runtime_error(ss.str());                                  \
    }                                                                      \
  }

class XFDTDCudaTensorHDException : public XFDTDException {
 public:
  XFDTDCudaTensorHDException(const std::string &mes) : XFDTDException{mes} {}
};

template <typename T, SizeType N>
class TensorHD {
 public:
 public:
  using DeviceTensor = Tensor<T, N>;

  template <typename WrappedTensor>
  TensorHD(const WrappedTensor &wrapped_tensor) {
    const auto &dim = wrapped_tensor.dimension();
    const auto &shape = wrapped_tensor.shape();
    const auto &stride = wrapped_tensor.strides();
    const auto size = wrapped_tensor.size();
    // remove const
    auto data = const_cast<T *>(wrapped_tensor.data());

    // check dim
    if (dim != DeviceTensor::dimension()) {
      std::stringstream ss;
      ss << "Wrong dim";
      throw XFDTDCudaTensorHDException{ss.str()};
    }

    for (SizeType i = 0; i < dim; ++i) {
      _shape[i] = shape[i];
    }

    _host_data = data;
  }

  // TensorHD(FixedArray<SizeType, N> shape) : _host{new HostTensor{shape}} {}

  ~TensorHD() {
    if (_device) {
      // don't call destructor to avoid double free
      // __destroyDeviceObject(_device);
      {
        auto err = cudaFree(_device);
        if (err != cudaSuccess) {
          std::cerr << "Failed to free device memory \n";
          std::abort();
        }
      }
      _device = nullptr;
    }

    if (_device_data) {
      {
        auto err = cudaFree(_device_data);
        if (err != cudaSuccess) {
          std::cerr << "Failed to free device memory \n";
          std::abort();
        }
      }
      _device_data = nullptr;
    }
  }

  XFDTD_CORE_CUDA_DUAL auto device() -> DeviceTensor * { return _device; }

  XFDTD_CORE_CUDA_DUAL auto device() const -> const DeviceTensor * {
    return _device;
  }

  auto hostData() { return _host_data; }

  auto hostData() const { return _host_data; }

  /**
   * @brief Reset host data. It doesn't free previous host data and copy data to
   * device. The only thing it does is to change the host data pointer.
   */
  auto resetHostData(T *data) { _host_data = data; }

  auto copyHostToDevice() -> void {
    if (_host_data == nullptr) {
      throw XFDTDCudaTensorHDException{"Host memory is not allocated"};
    }

    if (_device) {
      // Free previous device memory
      XFDTD_CORE_CUDA_CHECK_CUDA_ERROR(cudaFree(_device));
      _device = nullptr;
    }

    if (_device_data != nullptr) {
      XFDTD_CORE_CUDA_CHECK_CUDA_ERROR(cudaFree(_device_data));
      _device_data = nullptr;
    }

    // Malloc device
    XFDTD_CORE_CUDA_CHECK_CUDA_ERROR(
        cudaMalloc(&_device, sizeof(DeviceTensor)));

    // TODO: first: maclloc device data, second: malloc metadata
    // Copy tensor metadata
    auto host_tensor_matedata = DeviceTensor{};
    host_tensor_matedata._shape = _shape;
    host_tensor_matedata._strides = host_tensor_matedata.makeStride(_shape);
    host_tensor_matedata._size = host_tensor_matedata.makeSize(_shape);
    host_tensor_matedata._data = nullptr;
    XFDTD_CORE_CUDA_CHECK_CUDA_ERROR(cudaMemcpy(_device, &host_tensor_matedata,
                                                sizeof(DeviceTensor),
                                                cudaMemcpyHostToDevice));

    // Malloc device data
    XFDTD_CORE_CUDA_CHECK_CUDA_ERROR(
        cudaMalloc(&_device_data, host_tensor_matedata.size() * sizeof(T)));

    // Copy data
    XFDTD_CORE_CUDA_CHECK_CUDA_ERROR(cudaMemcpy(
        _device_data, _host_data, host_tensor_matedata.size() * sizeof(T),
        cudaMemcpyHostToDevice))

    // set decive data
    __kernelSetDeviceArrayData<<<1, 1>>>(_device, _device_data);
    cudaDeviceSynchronize();
  }

  auto copyDeviceToHost() -> void {
    if (!_device || !_device_data) {
      throw std::runtime_error("Device memory is not allocated");
    }

    // recieve meta data
    auto host_tensor_matedata = DeviceTensor{};
    XFDTD_CORE_CUDA_CHECK_CUDA_ERROR(cudaMemcpy(&host_tensor_matedata, _device,
                                                sizeof(DeviceTensor),
                                                cudaMemcpyDeviceToHost));
    host_tensor_matedata._data = nullptr;  // can't receive data in device
    // assume that shape will be never changed. Just copy
    XFDTD_CORE_CUDA_CHECK_CUDA_ERROR(
        cudaMemcpy(_host_data, _device_data,
                   host_tensor_matedata.size() * sizeof(T),
                   cudaMemcpyDeviceToHost););
  }

 protected:
 private:
  DeviceTensor *_device{};

  T *_device_data{};
  T *_host_data{};

  FixedArray<SizeType, N> _shape;
};

template <typename WrappedTensor, SizeType N>
class TensorHDWrapped : public TensorHD<typename WrappedTensor::value_type, N> {
 public:
 public:
  using T = typename WrappedTensor::value_type;
  TensorHDWrapped(WrappedTensor tensor)
      : TensorHD<T, N>{tensor}, _host_tensor{std::move(tensor)} {}

  auto tensor() const -> const WrappedTensor & { return _host_tensor; }

  auto tensor() -> WrappedTensor & { return _host_tensor; }

 private:
  WrappedTensor _host_tensor;
};

}  // namespace cuda

}  // namespace xfdtd

#endif  // __XFDTD_CORE_CUDA_TENSOR_HD_CUH__
