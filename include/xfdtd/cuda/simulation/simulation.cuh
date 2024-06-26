#ifndef __XFDTD_CORE_CUDA_SIMULATION_CUH__
#define __XFDTD_CORE_CUDA_SIMULATION_CUH__

#include <memory>

#include "xfdtd/cuda/common.cuh"
#include "xfdtd/cuda/grid_space/grid_space.cuh"
#include "xfdtd/grid_space/grid_space.h"
namespace xfdtd {

namespace cuda {

class Simulation {
 public:
  Simulation();

  auto init() -> void;

 private:
  std::shared_ptr<GridSpace> _grid_space;
  GridSpaceHD _grid_space_hd;
};

}  // namespace cuda

}  // namespace xfdtd

#endif  // __XFDTD_CORE_CUDA_SIMULATION_CUH__
