#include <xfdtd/monitor/monitor.h>

#include <filesystem>
#include <xtensor/xnpy.hpp>

namespace xfdtd {

Monitor::Monitor(std::unique_ptr<Shape> shape, std::string name,
                 std::string output_dir)
    : _shape{std::move(shape)},
      _name(std::move(name)),
      _output_dir(std::move(output_dir)) {}

const std::string& Monitor::name() const { return _name; }

const std::string& Monitor::outputDir() const { return _output_dir; }

const xt::xarray<double>& Monitor::data() const { return _data; }

xt::xarray<double>& Monitor::data() { return _data; }

void Monitor::setName(std::string name) { _name = std::move(name); }

void Monitor::setOutputDir(std::string output_dir) {
  _output_dir = std::move(output_dir);
}

void Monitor::output() {
  auto out_dir{std::filesystem::path(_output_dir)};
  if (!std::filesystem::exists(out_dir)) {
    std::filesystem::create_directories(out_dir);
  }

  auto out_file{out_dir / (name() + ".npy")};
  xt::dump_npy(out_file.string(), _data);
}

void Monitor::defaultInit(
    std::shared_ptr<const GridSpace> grid_space,
    std::shared_ptr<const CalculationParam> calculation_param,
    std::shared_ptr<const EMF> emf) {
  _grid_space = std::move(grid_space);
  _calculation_param = std::move(calculation_param);
  _emf = std::move(emf);
  if (_shape == nullptr) {
    return;
  }
  _grid_box = std::make_unique<GridBox>(_grid_space->getGridBox(_shape.get()));
}

const GridSpace* Monitor::gridSpacePtr() const { return _grid_space.get(); }

const CalculationParam* Monitor::calculationParamPtr() const {
  return _calculation_param.get();
}

const EMF* Monitor::emfPtr() const { return _emf.get(); }

const GridBox* Monitor::gridBoxPtr() const { return _grid_box.get(); }

}  // namespace xfdtd
