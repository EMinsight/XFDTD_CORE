#include <xfdtd/common/type_define.h>
#include <xfdtd/grid_space/grid_space.h>
#include <xfdtd/monitor/monitor.h>
#include <xfdtd/parallel/mpi_support.h>
#include <xfdtd/parallel/parallelized_config.h>

#include <algorithm>
#include <filesystem>
#include <xtensor/xio.hpp>
#include <xtensor/xnpy.hpp>

namespace xfdtd {

Monitor::Monitor(std::unique_ptr<Shape> shape, std::string name,
                 std::string output_dir)
    : _shape{std::move(shape)},
      _name(std::move(name)),
      _output_dir(std::move(output_dir)) {}

const std::unique_ptr<Shape>& Monitor::shape() const { return _shape; }

const std::string& Monitor::name() const { return _name; }

const std::string& Monitor::outputDir() const { return _output_dir; }

const Array<Real>& Monitor::data() const { return _data; }

Array<Real>& Monitor::data() { return _data; }

std::unique_ptr<Shape>& Monitor::shape() { return _shape; }

auto Monitor::globalTask() const -> IndexTask { return _global_task; }

auto Monitor::nodeTask() const -> IndexTask { return _node_task; }

void Monitor::setName(std::string name) { _name = std::move(name); }

void Monitor::setOutputDir(std::string output_dir) {
  _output_dir = std::move(output_dir);
}

void Monitor::output() {
  if (!valid() || !monitorMpiConfig().isRoot()) {
    return;
  }

  auto out_dir{std::filesystem::path(_output_dir)};
  if (!std::filesystem::exists(out_dir)) {
    std::filesystem::create_directories(out_dir);
  }

  auto out_file{out_dir / (name() + ".npy")};
  xt::dump_npy(out_file.string(), _data);
}

void Monitor::initTimeDependentVariable() {}

GridBox Monitor::globalGridBox() const { return _global_grid_box; }

GridBox Monitor::nodeGridBox() const { return _node_grid_box; }

const MpiConfig& Monitor::monitorMpiConfig() const {
  return _monitor_mpi_config;
}

std::string Monitor::toString() const {
  std::stringstream ss;
  ss << "Monitor: " << name() << "\n";
  ss << " Shape: " << shape()->toString() << "\n";
  ss << " Global : " << globalGridBox().toString() << "\n";
  ss << " Node : " << nodeGridBox().toString() << "\n";
  ss << " Node Grid Box in Global Grid Box: "
     << _grid_space->transformNodeToGlobal(nodeGridBox()).toString() << "\n";
  ss << " Output Dir: " << outputDir() << "\n";
  ss << " Shape of Data: " << xt::adapt(data().shape());
  ss << "\n " << monitorMpiConfig().toString() << "\n";
  return ss.str();
}

auto Monitor::makeMpiSubComm() -> void {
  auto& mpi_support = MpiSupport::instance();
  auto arr = std::vector<int>(mpi_support.size(), 0);

  int is_valid = static_cast<int>(valid());

  mpi_support.allGather(mpi_support.config(), &is_valid, sizeof(int),
                        arr.data(), sizeof(int));
#if defined(XFDTD_CORE_WITH_MPI)
  int color = (arr[mpi_support.rank()] == 1) ? (0) : (MPI_UNDEFINED);
#else
  int color = 0;
#endif

  int counter =
      std::count_if(arr.begin(), arr.end(), [](const auto& v) { return v; });

  _monitor_mpi_config =
      MpiConfig::makeSub(mpi_support.config(), color, counter);
}

auto Monitor::gatherData() -> void {}

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

  _global_grid_box = _grid_space->globalGridSpace()->getGridBox(shape().get());
  _node_grid_box = _grid_space->getGridBoxWithoutCheck(shape().get());

  _global_task = makeIndexTask(
      makeIndexRange(_global_grid_box.origin().i(), _global_grid_box.end().i()),
      makeIndexRange(_global_grid_box.origin().j(), _global_grid_box.end().j()),
      makeIndexRange(_global_grid_box.origin().k(),
                     _global_grid_box.end().k()));

  _node_task = makeIndexTask(
      makeIndexRange(_node_grid_box.origin().i(), _node_grid_box.end().i()),
      makeIndexRange(_node_grid_box.origin().j(), _node_grid_box.end().j()),
      makeIndexRange(_node_grid_box.origin().k(), _node_grid_box.end().k()));
}

auto Monitor::setGlobalTask(IndexTask task) -> void { _global_task = task; }

auto Monitor::setNodeTask(IndexTask task) -> void { _node_task = task; }

const GridSpace* Monitor::gridSpacePtr() const { return _grid_space.get(); }

const CalculationParam* Monitor::calculationParamPtr() const {
  return _calculation_param.get();
}

const EMF* Monitor::emfPtr() const { return _emf.get(); }

MpiConfig& Monitor::monitorMpiConfig() { return _monitor_mpi_config; }

auto Monitor::setGlobalGridBox(GridBox grid_box) -> void {
  _global_grid_box = grid_box;
}

auto Monitor::setNodeGridBox(GridBox grid_box) -> void {
  _node_grid_box = grid_box;
}

auto Monitor::nodeGridBox() -> GridBox& { return _node_grid_box; }

auto Monitor::nodeTask() -> IndexTask& { return _node_task; }

auto Monitor::valid() const -> bool { return nodeTask().valid(); }

}  // namespace xfdtd
