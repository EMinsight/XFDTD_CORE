#ifndef _XFDTD_CORE_CALCULATION_PARAM_H_
#define _XFDTD_CORE_CALCULATION_PARAM_H_

#include <xfdtd/calculation_param/fdtd_update_coefficient.h>
#include <xfdtd/calculation_param/material_param.h>
#include <xfdtd/calculation_param/time_param.h>
#include <xfdtd/common/type_define.h>
#include <xfdtd/coordinate_system/coordinate_system.h>
#include <xfdtd/electromagnetic_field/electromagnetic_field.h>
#include <xfdtd/exception/exception.h>
#include <xfdtd/util/transform.h>

#include <memory>

namespace xfdtd {

class GridSpace;

class XFDTDCalculationParamException : public XFDTDException {
 public:
  explicit XFDTDCalculationParamException(std::string message);

  ~XFDTDCalculationParamException() override = default;
};

class CalculationParam {
 public:
  CalculationParam();

  ~CalculationParam() = default;

  const std::unique_ptr<TimeParam>& timeParam() const;

  const std::unique_ptr<MaterialParam>& materialParam() const;

  const std::unique_ptr<FDTDUpdateCoefficient>& fdtdCoefficient() const;

  std::unique_ptr<TimeParam>& timeParam();

  std::unique_ptr<MaterialParam>& materialParam();

  std::unique_ptr<FDTDUpdateCoefficient>& fdtdCoefficient();

  void generateMaterialSpaceParam(const GridSpace* grid_space);

  void calculateCoefficient(const GridSpace* grid_space);

  void setMaterialParam(std::unique_ptr<MaterialParam> material_param);

  void setTimeParam(std::unique_ptr<TimeParam> time_param);

 private:
  std::unique_ptr<TimeParam> _time_param;
  std::unique_ptr<MaterialParam> _material_param;
  std::unique_ptr<FDTDUpdateCoefficient> _fdtd_coefficient;
};

}  // namespace xfdtd

#endif  // _XFDTD_CORE_CALCULATION_PARAM_H_
