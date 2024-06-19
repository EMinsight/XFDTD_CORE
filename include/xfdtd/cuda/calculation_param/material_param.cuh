#ifndef __XFDTD_CORE_CUDA_MATERIAL_PARAM_CUH__
#define __XFDTD_CORE_CUDA_MATERIAL_PARAM_CUH__

#include <xfdtd/calculation_param/calculation_param.h>

#include <memory>

#include "xfdtd/cuda/common.cuh"

namespace xfdtd {

namespace cuda {

class MaterialParam {
  friend class MaterialParamHD;
};

class MaterialParamHD {
  using HostMaterialParam = xfdtd::MaterialParam;

 public:
  MaterialParamHD() = default;

  XFDTD_CORE_CUDA_DUAL MaterialParamHD(
      std::shared_ptr<HostMaterialParam> host_material_param)
      : _host_material_param{host_material_param} {}

 private:
  std::shared_ptr<HostMaterialParam> _host_material_param{};
};

}  // namespace cuda

}  // namespace xfdtd

#endif  // __XFDTD_CORE_CUDA_MATERIAL_PARAM_CUH__
