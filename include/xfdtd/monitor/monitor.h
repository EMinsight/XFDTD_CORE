#ifndef _XFDTD_LIB_MONITOR_H_
#define _XFDTD_LIB_MONITOR_H_

#include <xfdtd/calculation_param/calculation_param.h>
#include <xfdtd/electromagnetic_field/electromagnetic_field.h>
#include <xfdtd/grid_space/grid_space.h>
#include <xfdtd/shape/shape.h>

#include <memory>
#include <string>
#include <xtensor/xarray.hpp>

#include "xfdtd/exception/exception.h"

namespace xfdtd {

class XFDTDMonitorException : public XFDTDException {
 public:
  explicit XFDTDMonitorException(
      std::string message = "XFDTD Monitor Exception")
      : XFDTDException(std::move(message)) {}
};

class Monitor {
 public:
  explicit Monitor(std::unique_ptr<Shape> shape, std::string name = "monitor",
                   std::string output_dir = "xfdtd_output");

  Monitor(const Monitor&) = delete;

  Monitor(Monitor&&) noexcept = default;

  Monitor& operator=(const Monitor&) = delete;

  Monitor& operator=(Monitor&&) noexcept = default;

  virtual ~Monitor() = default;

  virtual void init(std::shared_ptr<const GridSpace> grid_space,
                    std::shared_ptr<const CalculationParam> calculation_param,
                    std::shared_ptr<const EMF> emf) = 0;

  virtual void update() = 0;

  const std::unique_ptr<Shape>& shape() const;

  const std::string& name() const;

  const std::string& outputDir() const;

  const xt::xarray<double>& data() const;

  std::unique_ptr<Shape>& shape();

  xt::xarray<double>& data();

  void setName(std::string name);

  void setOutputDir(std::string output_dir);

  virtual void output();

 protected:
  void defaultInit(std::shared_ptr<const GridSpace> grid_space,
                   std::shared_ptr<const CalculationParam> calculation_param,
                   std::shared_ptr<const EMF> emf);

  const GridSpace* gridSpacePtr() const;

  const CalculationParam* calculationParamPtr() const;

  const EMF* emfPtr() const;

  const GridBox* gridBoxPtr() const;

 private:
  std::unique_ptr<Shape> _shape;
  std::string _name;
  std::string _output_dir;
  xt::xarray<double> _data;

  std::shared_ptr<const GridSpace> _grid_space;
  std::shared_ptr<const CalculationParam> _calculation_param;
  std::shared_ptr<const EMF> _emf;

  std::unique_ptr<GridBox> _grid_box;
};

}  // namespace xfdtd

#endif  // _XFDTD_LIB_MONITOR_H_
