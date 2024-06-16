#ifndef __XFDTD_CORE_CUDA_TENSOR_HD_CUH__
#define __XFDTD_CORE_CUDA_TENSOR_HD_CUH__

#include <xfdtd/cuda/common.cuh>
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

template <typename T, SizeType N, typename WrappedTensor>
class TensorHD {
 public:
 public:
  using HostTensor = Tensor<T, N>;
  using DeviceTensor = Tensor<T, N>;

  TensorHD(WrappedTensor &wrapped_tensor) {
    auto dim = wrapped_tensor.dimension();
    auto shape = wrapped_tensor.shape();
    auto stride = wrapped_tensor.strides();
    auto size = wrapped_tensor.size();
    auto data = wrapped_tensor.data();


  }

  // TensorHD(Array<SizeType, N> shape) : _host{new HostTensor{shape}} {}

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

    delete _host;
    _host = nullptr;
  }

  XFDTD_CORE_CUDA_DUAL auto device() -> DeviceTensor * { return _device; }

  auto host() -> HostTensor * { return _host; }

  XFDTD_CORE_CUDA_DUAL auto device() const -> const DeviceTensor * {
    return _device;
  }

  auto host() const -> const HostTensor * { return _host; }

  auto copyHostToDevice() -> void {
    if (!_host) {
      throw std::runtime_error("Host memory is not allocated");
    }
    if (_device) {
      // Free previous device memory
      XFDTD_CORE_CUDA_CHECK_CUDA_ERROR(cudaFree(_device));
      _device = nullptr;
    }

    if (_device_data != nullptr) {
      XFDTD_CORE_CUDA_CHECK_CUDA_ERROR(cudaFree(_device_data));
    }

    XFDTD_CORE_CUDA_CHECK_CUDA_ERROR(
        cudaMalloc(&_device, sizeof(DeviceTensor)));

    // {
    //   // Allocate device memory
    //   auto err = cudaMalloc(&_device, sizeof(DeviceTensor));
    //   if (err != cudaSuccess) {
    //     throw std::runtime_error("Failed to allocate device memory" +
    //                              std::string(cudaGetErrorString(err)));
    //   }
    // }

    {
      // Copy tensor metadata
      auto err = cudaMemcpy(_device, _host, sizeof(HostTensor),
                            cudaMemcpyHostToDevice);
      if (err != cudaSuccess) {
        throw std::runtime_error("Failed to copy to device memory" +
                                 std::string(cudaGetErrorString(err)));
      }
    }

    {
      // Malloc device data
      auto err = cudaMalloc(&_device_data, _host->size() * sizeof(T));
      if (err != cudaSuccess) {
        throw std::runtime_error("Failed to allocate device memory" +
                                 std::string(cudaGetErrorString(err)));
      }
    }

    {
      // Copy tensor data
      auto err = cudaMemcpy(_device_data, _host->_data,
                            _host->size() * sizeof(T), cudaMemcpyHostToDevice);
      if (err != cudaSuccess) {
        throw std::runtime_error("Failed to copy to device memory" +
                                 std::string(cudaGetErrorString(err)));
      }
    }

    __kernelSetDeviceArrayData<<<1, 1>>>(_device, _device_data);
  }

  auto copyDeviceToHost() -> void {
    if (!_device || !_device_data) {
      throw std::runtime_error("Device memory is not allocated");
    }

    {
      auto err = cudaMemcpy(_host, _device, sizeof(DeviceTensor),
                            cudaMemcpyDeviceToHost);
      if (err != cudaSuccess) {
        throw std::runtime_error("Failed to copy to host memory" +
                                 std::string(cudaGetErrorString(err)));
      }
    }

    assert(_device_data == _host->_data);
    _host->_data = new T[_host->size()];

    {
      auto err = cudaMemcpy(_host->_data, _device_data,
                            _host->size() * sizeof(T), cudaMemcpyDeviceToHost);
      if (err != cudaSuccess) {
        throw std::runtime_error("Failed to copy data to host memory" +
                                 std::string(cudaGetErrorString(err)));
      }
    }

    // auto new_host = new HostTensor();

    // {
    //   auto err = cudaMemcpy(new_host, _device, sizeof(DeviceTensor),
    //                         cudaMemcpyDeviceToHost);
    //   if (err != cudaSuccess) {
    //     throw std::runtime_error("Failed to copy to host memory" +
    //                              std::string(cudaGetErrorString(err)));
    //   }
    // }

    // assert(_device_data == new_host->_data);
    // new_host->_data = new T[new_host->size()];

    // {
    //   auto err =
    //       cudaMemcpy(new_host->_data, _device_data,
    //                  new_host->size() * sizeof(T), cudaMemcpyDeviceToHost);
    //   if (err != cudaSuccess) {
    //     throw std::runtime_error("Failed to copy data to host memory" +
    //                              std::string(cudaGetErrorString(err)));
    //   }
    // }

    // delete _host;
    // _host = new_host;
  }

 protected:
 private:
  DeviceTensor *_device{};
  HostTensor *_host{};
  T *_device_data{};
};

}  // namespace cuda

}  // namespace xfdtd

#endif  // __XFDTD_CORE_CUDA_TENSOR_HD_CUH__
