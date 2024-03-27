#include <xfdtd/divider/divider.h>
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

const xt::xarray<double>& Monitor::data() const { return _data; }

xt::xarray<double>& Monitor::data() { return _data; }

std::unique_ptr<Shape>& Monitor::shape() { return _shape; }

auto Monitor::globalTask() const -> Divider::IndexTask { return _global_task; }

auto Monitor::nodeTask() const -> Divider::IndexTask { return _node_task; }

void Monitor::setName(std::string name) { _name = std::move(name); }

void Monitor::setOutputDir(std::string output_dir) {
  _output_dir = std::move(output_dir);
}

void Monitor::output() {
  if (!nodeTask().valid() || !monitorMpiConfig().isRoot()) {
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

auto Monitor::initParallelizedConfig() -> void {
  auto& mpi_support = MpiSupport::instance();
  auto arr = std::vector<int>(mpi_support.size(), 0);

  int is_valid = static_cast<int>(nodeTask().valid());

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

  if (is_valid == 0) {
    return;
  }

#if defined(XFDTD_CORE_WITH_MPI)
  const auto node_box_origin_in_global =
      _grid_space->globalBox().origin() + _node_grid_box.origin();
  const auto g_origin = _global_grid_box.origin();

  auto nx = _node_task.xRange().size();
  auto ny = _node_task.yRange().size();
  auto nz = _node_task.zRange().size();
  auto stride_elem = _global_task.zRange().size();
  auto stride_vec = _global_task.zRange().size() * _global_task.yRange().size();
  auto disp = (node_box_origin_in_global.i() - g_origin.i()) * stride_vec +
              (node_box_origin_in_global.j() - g_origin.j()) * stride_elem +
              (node_box_origin_in_global.k() - g_origin.k());

  auto p = MpiSupport::Block::Profile{
      static_cast<int>(nx),          static_cast<int>(ny),
      static_cast<int>(nz),          static_cast<int>(stride_vec),
      static_cast<int>(stride_elem), static_cast<int>(disp)};

  if (_monitor_mpi_config.isRoot()) {
    _profiles.resize(counter);
  }

  auto profile_type = MpiSupport::TypeGuard{};
  MPI_Type_contiguous(sizeof(MpiSupport::Block::Profile), MPI_CHAR,
                      &profile_type._type);
  MPI_Type_commit(&profile_type._type);

  _block = MpiSupport::Block::make(p);

  mpi_support.gather(monitorMpiConfig(), &p, 1, profile_type, _profiles.data(),
                     1, profile_type, monitorMpiConfig().root());

  if (_monitor_mpi_config.isRoot()) {
    _blocks_mpi.reserve(_profiles.size());
    for (const auto& profile : _profiles) {
      _blocks_mpi.emplace_back(MpiSupport::Block::make(profile));
    }
  }
#endif
}

auto Monitor::gatherData() -> void {
  if (!nodeTask().valid() || monitorMpiConfig().size() == 1) {
    return;
  }

  auto& mpi_support = MpiSupport::instance();
  if (monitorMpiConfig().isRoot()) {
    xt::xarray<double> recv_buffer = xt::zeros<double>(
        {_global_task.xRange().size(), _global_task.yRange().size(),
         _global_task.zRange().size()});
    for (int i = 1; i < monitorMpiConfig().size(); ++i) {
      mpi_support.iRecv(monitorMpiConfig(), recv_buffer.data(), 1,
                        _blocks_mpi[i], i, 0);
    }

    auto index =
        mpi_support.iSendRecv(monitorMpiConfig(), data().data(), data().size(),
                              monitorMpiConfig().rank(), 0, recv_buffer.data(),
                              1, _block, monitorMpiConfig().rank(), 0);
    mpi_support.waitAll();
    if (index == -1) {
      return;
    }

    data() = std::move(recv_buffer);
  } else {
    mpi_support.send(monitorMpiConfig(), data().data(),
                     sizeof(double) * data().size(), monitorMpiConfig().root(),
                     0);
  }
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

  _global_grid_box = _grid_space->globalGridSpace()->getGridBox(shape().get());
  _node_grid_box = _grid_space->getGridBoxWithoutCheck(shape().get());

  _global_task = Divider::makeIndexTask(
      Divider::makeIndexRange(_global_grid_box.origin().i(),
                              _global_grid_box.end().i()),
      Divider::makeIndexRange(_global_grid_box.origin().j(),
                              _global_grid_box.end().j()),
      Divider::makeIndexRange(_global_grid_box.origin().k(),
                              _global_grid_box.end().k()));

  _node_task = Divider::makeIndexTask(
      Divider::makeIndexRange(_node_grid_box.origin().i(),
                              _node_grid_box.end().i()),
      Divider::makeIndexRange(_node_grid_box.origin().j(),
                              _node_grid_box.end().j()),
      Divider::makeIndexRange(_node_grid_box.origin().k(),
                              _node_grid_box.end().k()));
}

auto Monitor::setGlobalTask(Divider::IndexTask task) -> void {
  _global_task = task;
}

auto Monitor::setNodeTask(Divider::IndexTask task) -> void {
  _node_task = task;
}

const GridSpace* Monitor::gridSpacePtr() const { return _grid_space.get(); }

const CalculationParam* Monitor::calculationParamPtr() const {
  return _calculation_param.get();
}

const EMF* Monitor::emfPtr() const { return _emf.get(); }

MpiConfig& Monitor::monitorMpiConfig() { return _monitor_mpi_config; }

auto Monitor::nodeGridBox() -> GridBox& { return _node_grid_box; }

auto Monitor::nodeTask() -> Divider::IndexTask& { return _node_task; }

auto Monitor::mpiBlock() -> MpiSupport::Block& { return _block; }

auto Monitor::mpiBlock() const -> const MpiSupport::Block& { return _block; }

auto Monitor::mpiBlockArray() -> std::vector<MpiSupport::Block>& {
  return _blocks_mpi;
}

auto Monitor::mpiBlockArray() const -> const std::vector<MpiSupport::Block>& {
  return _blocks_mpi;
}

}  // namespace xfdtd
