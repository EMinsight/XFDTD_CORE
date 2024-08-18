#ifndef __XFDTD_COREs_ADE_UPDATOR_H__
#define __XFDTD_COREs_ADE_UPDATOR_H__

#include <xfdtd/coordinate_system/coordinate_system.h>
#include <xfdtd/material/ade_method/ade_method.h>

#include "updator/basic_updator.h"

namespace xfdtd {

class ADEUpdator : public BasicUpdator {
 public:
  ADEUpdator(std::shared_ptr<const GridSpace> grid_space,
             std::shared_ptr<const CalculationParam> calculation_param,
             std::shared_ptr<EMF> emf, IndexTask task,
             std::shared_ptr<ADEMethodStorage> ade_method_storage);

  ~ADEUpdator() override = default;

  auto& storage() const { return _storage; }

  auto& storage() { return _storage; }

 private:
  std::shared_ptr<ADEMethodStorage> _storage{};
};

}  // namespace xfdtd

#endif  // __XFDTD_COREs_ADE_UPDATOR_H__
