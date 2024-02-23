#ifndef _XFDTD_LIB_BASIC_UPDATOR_H_
#define _XFDTD_LIB_BASIC_UPDATOR_H_

#include "updator/updator.h"

namespace xfdtd {
class BasicUpdator : public Updator {
 public:
  BasicUpdator(std::shared_ptr<const GridSpace> grid_space,
               std::shared_ptr<const CalculationParam> calculation_param,
               std::shared_ptr<EMF> emf);

  ~BasicUpdator() override = default;

  void updateH() override;
};

class BasicUpdatorTEM : public BasicUpdator {
 public:
  BasicUpdatorTEM(std::shared_ptr<const GridSpace> grid_space,
                  std::shared_ptr<const CalculationParam> calculation_param,
                  std::shared_ptr<EMF> emf);

  ~BasicUpdatorTEM() override = default;

  void updateE() override;
};

class BasicUpdatorTE : public BasicUpdator {
 public:
  BasicUpdatorTE(std::shared_ptr<const GridSpace> grid_space,
                 std::shared_ptr<const CalculationParam> calculation_param,
                 std::shared_ptr<EMF> emf);

  ~BasicUpdatorTE() override = default;

  void updateE() override;
};

class BasicUpdator3D : public BasicUpdator {
 public:
  BasicUpdator3D(std::shared_ptr<const GridSpace> grid_space,
                 std::shared_ptr<const CalculationParam> calculation_param,
                 std::shared_ptr<EMF> emf);

  ~BasicUpdator3D() override = default;

  void updateE() override;
};

}  // namespace xfdtd

#endif  // _XFDTD_LIB_BASIC_UPDATOR_H_
