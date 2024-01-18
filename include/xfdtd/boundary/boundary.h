#ifndef _XFDTD_LIB_BOUNDARY_H_
#define _XFDTD_LIB_BOUNDARY_H_

#include <xfdtd/calculation_param/calculation_param.h>
#include <xfdtd/electromagnetic_field/electromagnetic_field.h>
#include <xfdtd/grid_space/grid_space.h>

namespace xfdtd {

class XFDTDBoundaryException : public XFDTDException {
 public:
  explicit XFDTDBoundaryException(std::string message = "XFDTD Boundary Exception")
      : XFDTDException(std::move(message)) {}
};

class Boundary {
 public:
  Boundary() = default;

  Boundary(const Boundary&) = default;

  Boundary(Boundary&&) noexcept = default;

  Boundary& operator=(const Boundary&) = default;

  Boundary& operator=(Boundary&&) noexcept = default;

  virtual ~Boundary() = default;

  virtual void init(std::shared_ptr<const GridSpace> grid_space,
                    std::shared_ptr<CalculationParam> calculation_param,
                    std::shared_ptr<EMF> emf) = 0;

  virtual void correctMaterialSpace() = 0;

  virtual void correctUpdateCoefficient() = 0;

  virtual void correctE() = 0;

  virtual void correctH() = 0;

 protected:
  void defaultInit(std::shared_ptr<const GridSpace> grid_space,
                   std::shared_ptr<CalculationParam> calculation_param,
                   std::shared_ptr<EMF> emf);

  const GridSpace* gridSpacePtr() const;

  CalculationParam* calculationParamPtr() const;

  EMF* emfPtr() const;

 private:
  std::shared_ptr<const GridSpace> _grid_space;
  std::shared_ptr<CalculationParam> _calculation_param;
  std::shared_ptr<EMF> _emf;
};

};  // namespace xfdtd

#endif  // _XFDTD_LIB_BOUNDARY_H_
