#ifndef _XFDTD_LIB_UPDATOR_H_
#define _XFDTD_LIB_UPDATOR_H_

#include <memory>

#include "xfdtd/calculation_param/calculation_param.h"
#include "xfdtd/electromagnetic_field/electromagnetic_field.h"
#include "xfdtd/grid_space/grid_space.h"

namespace xfdtd {

class Updator {
 public:
  Updator(std::shared_ptr<const GridSpace> grid_space,
          std::shared_ptr<const CalculationParam> calculation_param,
          std::shared_ptr<EMF> emf);

  Updator(const Updator&) = default;

  Updator(Updator&&) noexcept = default;

  Updator& operator=(const Updator&) = default;

  Updator& operator=(Updator&&) noexcept = default;

  virtual ~Updator() = default;

  virtual void updateE() = 0;

  virtual void updateH() = 0;

 protected:
  std::shared_ptr<const GridSpace> _grid_space;
  std::shared_ptr<const CalculationParam> _calculation_param;
  std::shared_ptr<EMF> _emf;
};

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

#endif  // _XFDTD_LIB_UPDATOR_H_
