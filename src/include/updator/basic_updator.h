#ifndef _XFDTD_CORE_BASIC_UPDATOR_H_
#define _XFDTD_CORE_BASIC_UPDATOR_H_

#include <xfdtd/divider/divider.h>

#include "updator/updator.h"

namespace xfdtd {
class BasicUpdator : public Updator {
 public:
  BasicUpdator(std::shared_ptr<const GridSpace> grid_space,
               std::shared_ptr<const CalculationParam> calculation_param,
               std::shared_ptr<EMF> emf, Divider::IndexTask task);

  ~BasicUpdator() override = default;

  void updateH() override;
};

class BasicUpdatorTEM : public BasicUpdator {
 public:
  BasicUpdatorTEM(std::shared_ptr<const GridSpace> grid_space,
                  std::shared_ptr<const CalculationParam> calculation_param,
                  std::shared_ptr<EMF> emf, Divider::IndexTask task);

  ~BasicUpdatorTEM() override = default;

  void updateE() override;

 protected:
  void updateEEdge() override { throw std::runtime_error("Not implemented"); }

  void updateHEdge() override { throw std::runtime_error("Not implemented"); }
};

class BasicUpdatorTE : public BasicUpdator {
 public:
  BasicUpdatorTE(std::shared_ptr<const GridSpace> grid_space,
                 std::shared_ptr<const CalculationParam> calculation_param,
                 std::shared_ptr<EMF> emf, Divider::IndexTask task);

  ~BasicUpdatorTE() override = default;

  void updateE() override;

  std::string toString() const override;

 protected:
  void updateEEdge() override;

  void updateHEdge() override{};
};

class BasicUpdator3D : public BasicUpdator {
 public:
  BasicUpdator3D(std::shared_ptr<const GridSpace> grid_space,
                 std::shared_ptr<const CalculationParam> calculation_param,
                 std::shared_ptr<EMF> emf, Divider::IndexTask task);

  ~BasicUpdator3D() override = default;

  std::string toString() const override;

  void updateE() override;

 protected:
  void updateEEdge() override;

  void updateHEdge() override {}

 private:
  // void updateExYN();

  // void updateExZN();

  // void updateExCornerYZ();

  // void updateEyZN();

  // void updateEyXN();

  // void updateEyCornerXZ();

  // void updateEzXN();

  // void updateEzYN();

  // void updateEzCornerXY();
};

}  // namespace xfdtd

#endif  // _XFDTD_CORE_BASIC_UPDATOR_H_
