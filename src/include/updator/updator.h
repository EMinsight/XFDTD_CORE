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

}  // namespace xfdtd

#endif  // _XFDTD_LIB_UPDATOR_H_
