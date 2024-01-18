#ifndef _XFDTD_LIB_OBJECT_H_
#define _XFDTD_LIB_OBJECT_H_

#include <xfdtd/grid_space/grid_space.h>

#include <memory>

#include "xfdtd/calculation_param/calculation_param.h"
#include "xfdtd/electromagnetic_field/electromagnetic_field.h"
#include "xfdtd/material/material.h"
#include "xfdtd/shape/shape.h"

namespace xfdtd {

class XFDTDObjectException : public XFDTDException {
 public:
  explicit XFDTDObjectException(std::string message = "XFDTD Object Exception")
      : XFDTDException(std::move(message)) {}
};

class Object {
 public:
  Object(std::string name, std::unique_ptr<Shape> shape,
         std::unique_ptr<Material> material);

  Object(const Object&) = delete;

  Object(Object&&) noexcept = default;

  Object& operator=(const Object&) = delete;

  Object& operator=(Object&&) noexcept = default;

  virtual ~Object() = default;

  virtual std::string toString() const;

  virtual void init(std::shared_ptr<const GridSpace> grid_space,
                    std::shared_ptr<CalculationParam> calculation_param,
                    std::shared_ptr<EMF> emf);

  virtual void correctMaterialSpace();

  virtual void correctUpdateCoefficient();

  virtual void correctE();

  virtual void correctH();

  std::string name() const;

  const std::unique_ptr<Shape>& shape() const;

  const std::unique_ptr<Material>& material() const;

 protected:
  void defaultCorrectMaterialSpace();

  Shape* shapePtr();

  Material* materialPtr();

  void handleDispersion();

  const GridSpace* gridSpacePtr() const;

  CalculationParam* calculationParamPtr();

  EMF* emfPtr();

  GridBox* gridBoxPtr() const;

 private:
  std::string _name;
  std::unique_ptr<Shape> _shape;
  std::unique_ptr<Material> _material;

  std::shared_ptr<const GridSpace> _grid_space;
  std::shared_ptr<CalculationParam> _calculation_param;
  std::shared_ptr<EMF> _emf;

  std::unique_ptr<GridBox> _grid_box;
};

}  // namespace xfdtd

#endif  // _XFDTD_LIB_OBJECT_H_
