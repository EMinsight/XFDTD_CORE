#ifndef __XFDTD_CORE_CUDA_FDTD_COEFFICIENT_CUH__
#define __XFDTD_CORE_CUDA_FDTD_COEFFICIENT_CUH__

#include <xfdtd/calculation_param/calculation_param.h>
#include <xfdtd/common/type_define.h>

#include <xfdtd/cuda/common.cuh>
#include <xfdtd/cuda/tensor.cuh>
#include <xfdtd/cuda/tensor_hd.cuh>

namespace xfdtd {

namespace cuda {

class FDTDCoefficient {
 public:
  friend class FDTDCoefficientHD;

 public:
  FDTDCoefficient() = default;

 public:
  XFDTD_CORE_CUDA_DUAL auto cexe() const -> const Array3D<Real>& {
    return *_cexe;
  }
  XFDTD_CORE_CUDA_DUAL auto cexhy() const -> const Array3D<Real>& {
    return *_cexhy;
  }
  XFDTD_CORE_CUDA_DUAL auto cexhz() const -> const Array3D<Real>& {
    return *_cexhz;
  }
  XFDTD_CORE_CUDA_DUAL auto ceye() const -> const Array3D<Real>& {
    return *_ceye;
  }
  XFDTD_CORE_CUDA_DUAL auto ceyhz() const -> const Array3D<Real>& {
    return *_ceyhz;
  }
  XFDTD_CORE_CUDA_DUAL auto ceyhx() const -> const Array3D<Real>& {
    return *_ceyhx;
  }
  XFDTD_CORE_CUDA_DUAL auto ceze() const -> const Array3D<Real>& {
    return *_ceze;
  }
  XFDTD_CORE_CUDA_DUAL auto cezhx() const -> const Array3D<Real>& {
    return *_cezhx;
  }
  XFDTD_CORE_CUDA_DUAL auto cezhy() const -> const Array3D<Real>& {
    return *_cezhy;
  }

  XFDTD_CORE_CUDA_DUAL auto chxh() const -> const Array3D<Real>& {
    return *_chxh;
  }
  XFDTD_CORE_CUDA_DUAL auto chxey() const -> const Array3D<Real>& {
    return *_chxey;
  }
  XFDTD_CORE_CUDA_DUAL auto chxez() const -> const Array3D<Real>& {
    return *_chxez;
  }
  XFDTD_CORE_CUDA_DUAL auto chyh() const -> const Array3D<Real>& {
    return *_chyh;
  }
  XFDTD_CORE_CUDA_DUAL auto chyez() const -> const Array3D<Real>& {
    return *_chyez;
  }
  XFDTD_CORE_CUDA_DUAL auto chyex() const -> const Array3D<Real>& {
    return *_chyex;
  }
  XFDTD_CORE_CUDA_DUAL auto chzh() const -> const Array3D<Real>& {
    return *_chzh;
  }
  XFDTD_CORE_CUDA_DUAL auto chzex() const -> const Array3D<Real>& {
    return *_chzex;
  }
  XFDTD_CORE_CUDA_DUAL auto chzey() const -> const Array3D<Real>& {
    return *_chzey;
  }

  XFDTD_CORE_CUDA_DUAL auto cexe() -> Array3D<Real>& { return *_cexe; }
  XFDTD_CORE_CUDA_DUAL auto cexhy() -> Array3D<Real>& { return *_cexhy; }
  XFDTD_CORE_CUDA_DUAL auto cexhz() -> Array3D<Real>& { return *_cexhz; }
  XFDTD_CORE_CUDA_DUAL auto ceye() -> Array3D<Real>& { return *_ceye; }
  XFDTD_CORE_CUDA_DUAL auto ceyhz() -> Array3D<Real>& { return *_ceyhz; }
  XFDTD_CORE_CUDA_DUAL auto ceyhx() -> Array3D<Real>& { return *_ceyhx; }
  XFDTD_CORE_CUDA_DUAL auto ceze() -> Array3D<Real>& { return *_ceze; }
  XFDTD_CORE_CUDA_DUAL auto cezhx() -> Array3D<Real>& { return *_cezhx; }
  XFDTD_CORE_CUDA_DUAL auto cezhy() -> Array3D<Real>& { return *_cezhy; }

  XFDTD_CORE_CUDA_DUAL auto chxh() -> Array3D<Real>& { return *_chxh; }
  XFDTD_CORE_CUDA_DUAL auto chxey() -> Array3D<Real>& { return *_chxey; }
  XFDTD_CORE_CUDA_DUAL auto chxez() -> Array3D<Real>& { return *_chxez; }
  XFDTD_CORE_CUDA_DUAL auto chyh() -> Array3D<Real>& { return *_chyh; }
  XFDTD_CORE_CUDA_DUAL auto chyez() -> Array3D<Real>& { return *_chyez; }
  XFDTD_CORE_CUDA_DUAL auto chyex() -> Array3D<Real>& { return *_chyex; }
  XFDTD_CORE_CUDA_DUAL auto chzh() -> Array3D<Real>& { return *_chzh; }
  XFDTD_CORE_CUDA_DUAL auto chzex() -> Array3D<Real>& { return *_chzex; }
  XFDTD_CORE_CUDA_DUAL auto chzey() -> Array3D<Real>& { return *_chzey; }

 private:
  Array3D<Real>* _cexe{};
  Array3D<Real>* _cexhy{};
  Array3D<Real>* _cexhz{};
  Array3D<Real>* _ceye{};
  Array3D<Real>* _ceyhz{};
  Array3D<Real>* _ceyhx{};
  Array3D<Real>* _ceze{};
  Array3D<Real>* _cezhx{};
  Array3D<Real>* _cezhy{};

  Array3D<Real>* _chxh{};
  Array3D<Real>* _chxey{};
  Array3D<Real>* _chxez{};
  Array3D<Real>* _chyh{};
  Array3D<Real>* _chyez{};
  Array3D<Real>* _chyex{};
  Array3D<Real>* _chzh{};
  Array3D<Real>* _chzex{};
  Array3D<Real>* _chzey{};
};

class FDTDCoefficientHD {
 public:
  using HostFDTDCoefficient = xfdtd::FDTDUpdateCoefficient;
  using DeviceFDTDCoefficient = xfdtd::cuda::FDTDCoefficient;

  FDTDCoefficientHD(HostFDTDCoefficient* host_fdtd_coefficient);

 public:
  auto copyHostToDevice() -> void;

  auto copyDeviceToHost() -> void;

  auto releaseDevice() -> void;

  auto host() { return _host_fdtd_coefficient; }

  auto device() { return _device_fdtd_coefficient; }

 private:
  HostFDTDCoefficient* _host_fdtd_coefficient{};
  DeviceFDTDCoefficient* _device_fdtd_coefficient{};

  TensorHD<Real, 3> _cexe_hd, _cexhy_hd, _cexhz_hd, _ceye_hd, _ceyhz_hd,
      _ceyhx_hd, _ceze_hd, _cezhx_hd, _cezhy_hd;
  TensorHD<Real, 3> _chxh_hd, _chxey_hd, _chxez_hd, _chyh_hd, _chyez_hd,
      _chyex_hd, _chzh_hd, _chzex_hd, _chzey_hd;
};

}  // namespace cuda

}  // namespace xfdtd

#endif  // __XFDTD_CORE_CUDA_FDTD_COEFFICIENT_CUH__
