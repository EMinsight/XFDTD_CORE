#ifndef _XFDTD_CORE_OBJECT_H_
#define _XFDTD_CORE_OBJECT_H_

#include <xfdtd/calculation_param/calculation_param.h>
#include <xfdtd/common/index_task.h>
#include <xfdtd/common/type_define.h>
#include <xfdtd/electromagnetic_field/electromagnetic_field.h>
#include <xfdtd/grid_space/grid_space.h>
#include <xfdtd/material/material.h>
#include <xfdtd/shape/shape.h>

#include <memory>

namespace xfdtd {

class ADEMethodStorage;

class XFDTDObjectException : public XFDTDException {
 public:
  explicit XFDTDObjectException(std::string message = "XFDTD Object Exception")
      : XFDTDException(std::move(message)) {}
};

class Corrector;

class Object {
 public:
  Object(std::string name, std::unique_ptr<Shape> shape,
         std::shared_ptr<Material> material);

  Object(const Object&) = delete;

  Object(Object&&) noexcept = default;

  Object& operator=(const Object&) = delete;

  Object& operator=(Object&&) noexcept = default;

  virtual ~Object();

  virtual std::string toString() const;

  virtual void init(std::shared_ptr<const GridSpace> grid_space,
                    std::shared_ptr<CalculationParam> calculation_param,
                    std::shared_ptr<EMF> emf);

  virtual void correctMaterialSpace(Index index);

  virtual void correctUpdateCoefficient();

  virtual void handleDispersion(
      std::shared_ptr<ADEMethodStorage> ade_method_storage);

  virtual void initTimeDependentVariable();

  virtual void correctE();

  virtual void correctH();

  virtual std::unique_ptr<Corrector> generateCorrector(const Task<Index>& task);

  std::string name() const;

  const std::unique_ptr<Shape>& shape() const;

  const auto& material() const;

  auto&& material();

  const std::shared_ptr<const GridSpace>& gridSpace() const {
    return _grid_space;
  }

  std::shared_ptr<CalculationParam> calculationParam() const {
    return _calculation_param;
  }

  std::shared_ptr<EMF> emf() const { return _emf; }

  auto materialIndex() const -> Index;

 protected:

  auto setMaterialIndex(Index index) -> void;

  void defaultCorrectMaterialSpace(Index index = -1);

  Shape* shapePtr();

  Material* materialPtr();

  const GridSpace* gridSpacePtr() const;

  CalculationParam* calculationParamPtr();

  EMF* emfPtr();

  GridBox* gridBoxPtr() const;

  GridBox globalGridBox() const;

 private:
  std::string _name;
  std::unique_ptr<Shape> _shape;
  std::shared_ptr<Material> _material;

  std::shared_ptr<const GridSpace> _grid_space;
  std::shared_ptr<CalculationParam> _calculation_param;
  std::shared_ptr<EMF> _emf;

  Index _material_index;
  std::unique_ptr<GridBox> _grid_box;
  GridBox _global_grid_box;
};

inline const auto& Object::material() const { return _material; }

inline auto&& Object::material() { return _material; }

}  // namespace xfdtd

#endif  // _XFDTD_CORE_OBJECT_H_
