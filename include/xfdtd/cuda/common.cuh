#ifndef __XFDTD_CORE_CUDA_COMMON_CUH__
#define __XFDTD_CORE_CUDA_COMMON_CUH__

#include <cstddef>

namespace xfdtd::cuda {

#ifdef XFDTD_CORE_WITH_CUDA
#define XFDTD_CORE_CUDA_GLOBAL __global__
#define XFDTD_CORE_CUDA_HOST __host__
#define XFDTD_CORE_CUDA_DEVICE __device__
#define XFDTD_CORE_CUDA_DUAL __host__ __device__
#else
#define XFDTD_CORE_CUDA_GLOBAL
#define XFDTD_CORE_CUDA_HOST
#define XFDTD_CORE_CUDA_DEVICE
#define XFDTD_CORE_CUDA_DUAL
#endif

using SizeType = std::size_t;

// forward declaration
template <typename T, SizeType N>
class Tensor;

template <typename T>
using Array1D = Tensor<T, 1>;

template <typename T>
using Array2D = Tensor<T, 2>;

template <typename T>
using Array3D = Tensor<T, 3>;

template <typename T>
using Array4D = Tensor<T, 4>;

}  // namespace xfdtd::cuda

#endif  // __XFDTD_CORE_CUDA_COMMON_CUH__
