#include <xfdtd/common/constant.h>
#include <xfdtd/electromagnetic_field/electromagnetic_field.h>
#include <xfdtd/grid_space/grid_space.h>
#include <xfdtd/nffft/nffft.h>
#include <xfdtd/parallel/mpi_config.h>
#include <xfdtd/parallel/mpi_support.h>
#include <xfdtd/util/transform.h>

#include <memory>
#include <tuple>
#include <utility>
#include <vector>
#include <xtensor.hpp>
#include <xtensor/xcomplex.hpp>
#include <xtensor/xnpy.hpp>

namespace xfdtd {

NFFFT::NFFFT(Index distance_x, Index distance_y, Index distance_z)
    : _distance_x{distance_x},
      _distance_y{distance_y},
      _distance_z{distance_z} {
  if (distance_x == 0) {
    throw XFDTDNFFFTException("distance_x == 0");
  }

  if (distance_y == 0) {
    throw XFDTDNFFFTException("distance_y == 0");
  }

  if (distance_z == 0) {
    throw XFDTDNFFFTException("distance_z == 0");
  }
}

NFFFT::NFFFT(Index distance_x, Index distance_y, Index distance_z,
             std::string_view output_dir)
    : NFFFT(distance_x, distance_y, distance_z) {
  setOutputDir(output_dir);
}

NFFFT::~NFFFT() = default;

auto NFFFT::valid() const -> bool {
  return _node_task_surface_xn.valid() || _node_task_surface_xp.valid() ||
         _node_task_surface_yn.valid() || _node_task_surface_yp.valid() ||
         _node_task_surface_zn.valid() || _node_task_surface_zp.valid();
}

void NFFFT::defaultInit(
    std::shared_ptr<const GridSpace> grid_space,
    std::shared_ptr<const CalculationParam> calculation_param,
    std::shared_ptr<const EMF> emf) {
  _grid_space = std::move(grid_space);
  _calculation_param = std::move(calculation_param);
  _emf = std::move(emf);

  initGlobal();

  initNode();
}

auto NFFFT::initParallelizedConfig() -> void {
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

  _nffft_mpi_config = MpiConfig::makeSub(mpi_support.config(), color, counter);

  if (!valid() && mpi_support.size() <= 1) {
    throw XFDTDNFFFTException("Invalid task");
  }
}

auto NFFFT::initGlobal() -> void {
  // first, make global box
  const auto grid_space = gridSpace();
  const auto global_nx = grid_space->globalGridSpace()->sizeX();
  const auto global_ny = grid_space->globalGridSpace()->sizeY();
  const auto global_nz = grid_space->globalGridSpace()->sizeZ();

  const auto distance_x = distanceX();
  const auto distance_y = distanceY();
  const auto distance_z = distanceZ();

  const auto global_is = distance_x;
  const auto global_js = distance_y;
  const auto global_ks = distance_z;

  const auto global_ie = global_nx - distance_x;
  const auto global_je = global_ny - distance_y;
  const auto global_ke = global_nz - distance_z;

  if (global_ie <= global_is || global_je <= global_js ||
      global_ke <= global_ks) {
    throw XFDTDNFFFTException("Global box is not legal");
  }

  setGlobalBox(GridBox{Grid{global_is, global_js, global_ks},
                       Grid{global_ie - global_is, global_je - global_js,
                            global_ke - global_ks}});

  setGlobalTaskSurfaceXN(makeIndexTask(makeIndexRange(global_is, global_is + 1),
                                       makeIndexRange(global_js, global_je),
                                       makeIndexRange(global_ks, global_ke)));

  setGlobalTaskSurfaceXP(makeIndexTask(makeIndexRange(global_ie, global_ie + 1),
                                       makeIndexRange(global_js, global_je),
                                       makeIndexRange(global_ks, global_ke)));

  setGlobalTaskSurfaceYN(makeIndexTask(makeIndexRange(global_is, global_ie),
                                       makeIndexRange(global_js, global_js + 1),
                                       makeIndexRange(global_ks, global_ke)));

  setGlobalTaskSurfaceYP(makeIndexTask(makeIndexRange(global_is, global_ie),
                                       makeIndexRange(global_je, global_je + 1),
                                       makeIndexRange(global_ks, global_ke)));

  setGlobalTaskSurfaceZN(
      makeIndexTask(makeIndexRange(global_is, global_ie),
                    makeIndexRange(global_js, global_je),
                    makeIndexRange(global_ks, global_ks + 1)));

  setGlobalTaskSurfaceZP(
      makeIndexTask(makeIndexRange(global_is, global_ie),
                    makeIndexRange(global_js, global_je),
                    makeIndexRange(global_ke, global_ke + 1)));
}

auto NFFFT::initNode() -> void {
  const auto grid_space = gridSpace();
  auto cube = std::make_unique<Cube>(
      grid_space->globalGridSpace()->getGridOriginVector(globalBox().origin()),
      grid_space->globalGridSpace()->getGridOriginVector(globalBox().end()) -
          grid_space->globalGridSpace()->getGridOriginVector(
              globalBox().origin()));
  setNodeBox(grid_space->getGridBoxWithoutCheck(cube.get()));

  setNodeTaskSurfaceXN(makeNodeAxisTask(Axis::Direction::XN));
  setNodeTaskSurfaceXP(makeNodeAxisTask(Axis::Direction::XP));
  setNodeTaskSurfaceYN(makeNodeAxisTask(Axis::Direction::YN));
  setNodeTaskSurfaceYP(makeNodeAxisTask(Axis::Direction::YP));
  setNodeTaskSurfaceZN(makeNodeAxisTask(Axis::Direction::ZN));
  setNodeTaskSurfaceZP(makeNodeAxisTask(Axis::Direction::ZP));
}

auto NFFFT::makeNodeAxisTask(const Axis::Direction& direction) -> IndexTask {
  const auto global_box = globalBox();
  const auto node_box = nodeBox();
  const auto node_global_offset = gridSpace()->globalBox().origin();
  const auto node_lower = gridSpace()->box().origin();
  const auto node_upper = gridSpace()->box().end();

  const auto xyz = Axis::fromDirectionToXYZ(direction);

  const auto [global_offset_a, global_offset_b, global_offset_c] =
      transform::xYZToABC(
          std::tuple{node_global_offset.i(), node_global_offset.j(),
                     node_global_offset.k()},
          xyz);

  const auto [global_range_a, global_range_b, global_range_c] =
      transform::xYZToABC(
          std::tuple{
              makeIndexRange(global_box.origin().i(), global_box.end().i()),
              makeIndexRange(global_box.origin().j(), global_box.end().j()),
              makeIndexRange(global_box.origin().k(), global_box.end().k())},
          xyz);

  const auto [node_lower_a, node_lower_b, node_lower_c] = transform::xYZToABC(
      std::tuple{node_lower.i(), node_lower.j(), node_lower.k()}, xyz);

  const auto [node_upper_a, node_upper_b, node_upper_c] = transform::xYZToABC(
      std::tuple{node_upper.i(), node_upper.j(), node_upper.k()}, xyz);

  auto [node_range_a, node_range_b, node_range_c] = transform::xYZToABC(
      std::tuple{makeIndexRange(node_box.origin().i(), node_box.end().i()),
                 makeIndexRange(node_box.origin().j(), node_box.end().j()),
                 makeIndexRange(node_box.origin().k(), node_box.end().k())},
      xyz);

  auto invalid_task = IndexTask{};
  auto range_c = IndexRange{};

  const auto global_dest_c =
      (Axis::directionNegative(direction) ? (global_range_c.start())
                                          : (global_range_c.end()));
  if (node_range_c.start() + global_offset_c == global_dest_c) {
    if (node_range_c.start() <= node_lower_c) {
      // invalid
      return invalid_task;
    }

    if (node_range_c.start() == node_lower_c + 1) {
      // prefer to c_s = node_upper_c - 2
      return invalid_task;
    }

    range_c = makeIndexRange(node_range_c.start(), node_range_c.start() + 1);
  }

  if (node_range_c.end() + global_offset_c == global_dest_c) {
    if (node_range_c.end() >= node_upper_c) {
      // invalid
      return invalid_task;
    }

    range_c = makeIndexRange(node_range_c.end(), node_range_c.end() + 1);
  }

  if (!range_c.valid()) {
    return invalid_task;
  }

  auto make_secondary_axis_range = [](const auto& lower, const auto& upper,
                                      auto start, auto end) {
    if (start <= lower && lower == 0) {
      start = lower + 1;
    }

    if (end >= upper) {
      end = upper - 1;
    }

    return makeIndexRange(start, end);
  };

  auto range_a = make_secondary_axis_range(
      node_lower_a, node_upper_a, node_range_a.start(), node_range_a.end());
  auto range_b = make_secondary_axis_range(
      node_lower_b, node_upper_b, node_range_b.start(), node_range_b.end());

  auto [range_x, range_y, range_z] =
      transform::aBCToXYZ(std::tuple{range_a, range_b, range_c}, xyz);

  return makeIndexTask(range_x, range_y, range_z);
}

auto NFFFT::distanceX() const -> Index { return _distance_x; }

auto NFFFT::distanceY() const -> Index { return _distance_y; }

auto NFFFT::distanceZ() const -> Index { return _distance_z; }

auto NFFFT::outputDir() const -> std::string { return _output_dir; }

auto NFFFT::globalBox() const -> GridBox { return _global_box; }

auto NFFFT::nodeBox() const -> GridBox { return _node_box; }

auto NFFFT::globalTaskSurfaceXN() const -> IndexTask {
  return _global_task_surface_xn;
}

auto NFFFT::globalTaskSurfaceXP() const -> IndexTask {
  return _global_task_surface_xp;
}

auto NFFFT::globalTaskSurfaceYN() const -> IndexTask {
  return _global_task_surface_yn;
}

auto NFFFT::globalTaskSurfaceYP() const -> IndexTask {
  return _global_task_surface_yp;
}

auto NFFFT::globalTaskSurfaceZN() const -> IndexTask {
  return _global_task_surface_zn;
}

auto NFFFT::globalTaskSurfaceZP() const -> IndexTask {
  return _global_task_surface_zp;
}

auto NFFFT::nodeTaskSurfaceXN() const -> IndexTask {
  return _node_task_surface_xn;
}

auto NFFFT::nodeTaskSurfaceXP() const -> IndexTask {
  return _node_task_surface_xp;
}

auto NFFFT::nodeTaskSurfaceYN() const -> IndexTask {
  return _node_task_surface_yn;
}

auto NFFFT::nodeTaskSurfaceYP() const -> IndexTask {
  return _node_task_surface_yp;
}

auto NFFFT::nodeTaskSurfaceZN() const -> IndexTask {
  return _node_task_surface_zn;
}

auto NFFFT::nodeTaskSurfaceZP() const -> IndexTask {
  return _node_task_surface_zp;
}

auto NFFFT::setOutputDir(std::string_view output_dir) -> void {
  _output_dir = output_dir;
}

auto NFFFT::nffftMPIConfig() const -> const MpiConfig& {
  return _nffft_mpi_config;
}

auto NFFFT::setGlobalBox(GridBox box) -> void { _global_box = box; }

auto NFFFT::setNodeBox(GridBox box) -> void { _node_box = box; }

auto NFFFT::setGlobalTaskSurfaceXN(IndexTask task) -> void {
  _global_task_surface_xn = task;
}

auto NFFFT::setGlobalTaskSurfaceXP(IndexTask task) -> void {
  _global_task_surface_xp = task;
}

auto NFFFT::setGlobalTaskSurfaceYN(IndexTask task) -> void {
  _global_task_surface_yn = task;
}

auto NFFFT::setGlobalTaskSurfaceYP(IndexTask task) -> void {
  _global_task_surface_yp = task;
}

auto NFFFT::setGlobalTaskSurfaceZN(IndexTask task) -> void {
  _global_task_surface_zn = task;
}

auto NFFFT::setGlobalTaskSurfaceZP(IndexTask task) -> void {
  _global_task_surface_zp = task;
}

auto NFFFT::setNodeTaskSurfaceXN(IndexTask task) -> void {
  _node_task_surface_xn = task;
}

auto NFFFT::setNodeTaskSurfaceXP(IndexTask task) -> void {
  _node_task_surface_xp = task;
}

auto NFFFT::setNodeTaskSurfaceYN(IndexTask task) -> void {
  _node_task_surface_yn = task;
}

auto NFFFT::setNodeTaskSurfaceYP(IndexTask task) -> void {
  _node_task_surface_yp = task;
}

auto NFFFT::setNodeTaskSurfaceZN(IndexTask task) -> void {
  _node_task_surface_zn = task;
}

auto NFFFT::setNodeTaskSurfaceZP(IndexTask task) -> void {
  _node_task_surface_zp = task;
}

auto NFFFT::gridSpace() const -> std::shared_ptr<const GridSpace> {
  return _grid_space;
}

auto NFFFT::calculationParam() const
    -> std::shared_ptr<const CalculationParam> {
  return _calculation_param;
}

auto NFFFT::emf() const -> std::shared_ptr<const EMF> { return _emf; }

}  // namespace xfdtd
