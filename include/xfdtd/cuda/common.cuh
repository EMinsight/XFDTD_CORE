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

}  // namespace xfdtd::cuda

#endif  // __XFDTD_CORE_CUDA_COMMON_CUH__
