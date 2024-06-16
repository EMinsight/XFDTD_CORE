#ifndef __XFDTD_CORD_CUDA_MEMORY_CUH__
#define __XFDTD_CORD_CUDA_MEMORY_CUH__

#include <xfdtd/cuda/common.cuh>

namespace xfdtd {

namespace cuda {

template <typename D, typename T>
XFDTD_CORE_CUDA_GLOBAL auto __kernelSetDeviceArrayData(D *device,
                                                       T *data) -> void {
  device->setData(data);
}

}  // namespace cuda

}  // namespace xfdtd

#endif  // __XFDTD_CORD_CUDA_MEMORY_CUH__
