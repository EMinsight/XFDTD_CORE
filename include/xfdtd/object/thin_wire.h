#ifndef _XFDTD_LIB_THIN_WIRE_H_
#define _XFDTD_LIB_THIN_WIRE_H_

#include "xfdtd/coordinate_system/coordinate_system.h"
#include "xfdtd/object/object.h"
#include "xfdtd/shape/cylinder.h"

namespace xfdtd {

class ThinWire : public Object {
 public:
  ThinWire(const std::string& name, std::unique_ptr<Cylinder> shape);

  ThinWire(const ThinWire&) = delete;

  ThinWire(ThinWire&&) noexcept = default;

  ThinWire& operator=(const ThinWire&) = delete;

  ThinWire& operator=(ThinWire&&) noexcept = default;

  ~ThinWire() override = default;

  void init(std::shared_ptr<const GridSpace> grid_space,
            std::shared_ptr<CalculationParam> calculation_param,
            std::shared_ptr<EMF> emf) override;

  void correctUpdateCoefficient() override;

 private:
  Axis::XYZ _axis;
  double _radius;

  bool isEnoughThin() const;
};

}  // namespace xfdtd

#endif  // _XFDTD_LIB_THIN_WIRE_H_
