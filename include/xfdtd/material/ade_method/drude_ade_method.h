#ifndef __XFDTD_CORE_DRUDE_ADE_METHOD_H__
#define __XFDTD_CORE_DRUDE_ADE_METHOD_H__

#include <xfdtd/material/ade_method/ade_method.h>

namespace xfdtd {

class DrudeADEMethodStorage : public ADEMethodStorage {
 public:
  DrudeADEMethodStorage(Index num_pole, Index nx, Index ny, Index nz);

  auto correctCoeff(Index i, Index j, Index k,
                    const LinearDispersiveMaterial& linear_dispersive_material,
                    const std::shared_ptr<const GridSpace>& grid_space,
                    const std::shared_ptr<CalculationParam>& calculation_param)
      -> void override;

 private:
};

}  // namespace xfdtd

#endif  // __XFDTD_CORE_DRUDE_ADE_METHOD_H__
