#ifndef __XFDTD_CORE_M_LOR_ADE_METHOD_H__
#define __XFDTD_CORE_M_LOR_ADE_METHOD_H__

#include <xfdtd/material/ade_method/ade_method.h>

namespace xfdtd {

class MLorentzADEMethodStorage : public ADEMethodStorage {
 public:
  MLorentzADEMethodStorage(Index num_pole, Index nx, Index ny, Index nz);

  auto correctCoeff(Index i, Index j, Index k,
                    const LinearDispersiveMaterial& linear_dispersive_material,
                    const std::shared_ptr<const GridSpace>& grid_space,
                    const std::shared_ptr<CalculationParam>& calculation_param)
      -> void override;

 private:
};

}  // namespace xfdtd

#endif  // __XFDTD_CORE_M_LOR_ADE_METHOD_H__
