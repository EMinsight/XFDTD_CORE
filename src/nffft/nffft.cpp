#include <xfdtd/electromagnetic_field/electromagnetic_field.h>
#include <xfdtd/grid_space/grid_space.h>
#include <xfdtd/nffft/nffft.h>
#include <xfdtd/parallel/mpi_support.h>
#include <xfdtd/util/constant.h>
#include <xfdtd/util/transform.h>

#include <complex>
#include <cstdlib>
#include <filesystem>
#include <memory>
#include <sstream>
#include <tuple>
#include <utility>
#include <vector>
#include <xtensor.hpp>
#include <xtensor/xcomplex.hpp>
#include <xtensor/xnpy.hpp>

#include "nffft/nffft_fd_data.h"
#include "xfdtd/parallel/mpi_config.h"

namespace xfdtd {

NFFFT::NFFFT(std::size_t distance_x, std::size_t distance_y,
             std::size_t distance_z, xt::xarray<double> frequencies,
             std::string output_dir)
    : _distance_x{distance_x},
      _distance_y{distance_y},
      _distance_z{distance_z},
      _frequencies{std::move(frequencies)},
      _output_dir{std::move(output_dir)} {
  if (distance_x == 0) {
    throw XFDTDNFFFTException("distance_x == 0");
  }

  if (distance_y == 0) {
    throw XFDTDNFFFTException("distance_y == 0");
  }

  if (distance_z == 0) {
    throw XFDTDNFFFTException("distance_z == 0");
  }

  if (_frequencies.size() == 0) {
    throw XFDTDNFFFTException("frequencies.size() == 0");
  }
}

NFFFT::~NFFFT() = default;

auto NFFFT::valid() const -> bool {
  return _node_task_surface_xn.valid() || _node_task_surface_xp.valid() ||
         _node_task_surface_yn.valid() || _node_task_surface_yp.valid() ||
         _node_task_surface_zn.valid() || _node_task_surface_zp.valid();
}

void NFFFT::init(std::shared_ptr<const GridSpace> grid_space,
                 std::shared_ptr<const CalculationParam> calculation_param,
                 std::shared_ptr<const EMF> emf) {
  _grid_space = std::move(grid_space);
  _calculation_param = std::move(calculation_param);
  _emf = std::move(emf);

  initGlobal();

  initNode();

  if (!valid()) {
    return;
  }

  generateSurface();
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

void NFFFT::initTimeDependentVariable() {
  const auto total_time_step{
      calculationParamPtr()->timeParam()->endTimeStep() -
      calculationParamPtr()->timeParam()->startTimeStep()};
  const auto dt{calculationParamPtr()->timeParam()->dt()};
  for (auto&& fd_data : _fd_plane_data) {
    fd_data.initDFT(total_time_step, dt);
  }
}

auto NFFFT::update() -> void {
  auto current_time_step{calculationParamPtr()->timeParam()->currentTimeStep()};
  for (auto f{0}; f < _frequencies.size(); ++f) {
    for (auto&& fd_data : _fd_plane_data) {
      fd_data.update(current_time_step);
    }
  }
}

auto NFFFT::processFarField(const xt::xtensor<double, 1>& theta, double phi,
                            const std::string& sub_dir,
                            const Vector& origin) const -> void {
  processFarField(theta, xt::xtensor<double, 1>{phi}, sub_dir, origin);
}

auto NFFFT::processFarField(double theta, const xt::xtensor<double, 1>& phi,
                            const std::string& sub_dir,
                            const Vector& origin) const -> void {
  processFarField(xt::xtensor<double, 1>{theta}, phi, sub_dir, origin);
}

void NFFFT::outputRadiationPower() {
  if (!valid()) {
    return;
  }

  auto num_freq = _fd_plane_data.size();
  xt::xtensor<double, 1> freq_arr = xt::zeros<double>({num_freq});
  xt::xtensor<double, 1> node_power_arr = xt::zeros<double>({num_freq});
  xt::xtensor<double, 1> power_arr = xt::zeros<double>({num_freq});

  for (auto i = 0; i < num_freq; ++i) {
    freq_arr(i) = _fd_plane_data[i].frequency();
    node_power_arr(i) = _fd_plane_data[i].power();
  }

  MpiSupport::instance().reduceSum(nffftMPIConfig(), node_power_arr.data(),
                                   power_arr.data(), node_power_arr.size());

  if (!nffftMPIConfig().isRoot()) {
    return;
  }

  if (!std::filesystem::exists(_output_dir)) {
    std::filesystem::create_directories(_output_dir);
  }
  std::cout << "Rank: " << nffftMPIConfig().rank() << " do check power"
            << "\n";

  auto data = xt::stack(xt::xtuple(freq_arr, power_arr));
  xt::dump_npy((std::filesystem::path{_output_dir} / "power.npy").string(),
               data);
  std::cout << "Rank: " << nffftMPIConfig().rank() << " do write power done"
            << "\n";
}

auto NFFFT::initGlobal() -> void {
  // first, make global box
  const auto grid_space = gridSpacePtr();
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

  setGlobalTaskSurfaceXN(
      Divider::makeIndexTask(Divider::makeIndexRange(global_is, global_is + 1),
                             Divider::makeIndexRange(global_js, global_je),
                             Divider::makeIndexRange(global_ks, global_ke)));

  setGlobalTaskSurfaceXP(
      Divider::makeIndexTask(Divider::makeIndexRange(global_ie, global_ie + 1),
                             Divider::makeIndexRange(global_js, global_je),
                             Divider::makeIndexRange(global_ks, global_ke)));

  setGlobalTaskSurfaceYN(
      Divider::makeIndexTask(Divider::makeIndexRange(global_is, global_ie),
                             Divider::makeIndexRange(global_js, global_js + 1),
                             Divider::makeIndexRange(global_ks, global_ke)));

  setGlobalTaskSurfaceYP(
      Divider::makeIndexTask(Divider::makeIndexRange(global_is, global_ie),
                             Divider::makeIndexRange(global_je, global_je + 1),
                             Divider::makeIndexRange(global_ks, global_ke)));

  setGlobalTaskSurfaceZN(Divider::makeIndexTask(
      Divider::makeIndexRange(global_is, global_ie),
      Divider::makeIndexRange(global_js, global_je),
      Divider::makeIndexRange(global_ks, global_ks + 1)));

  setGlobalTaskSurfaceZP(Divider::makeIndexTask(
      Divider::makeIndexRange(global_is, global_ie),
      Divider::makeIndexRange(global_js, global_je),
      Divider::makeIndexRange(global_ke, global_ke + 1)));

  const auto& global_grid_space = grid_space->globalGridSpace();
  _cube =
      std::make_unique<Cube>(global_grid_space->getGridOriginVector(
                                 Grid{global_is, global_js, global_ks}),
                             global_grid_space->getGridOriginVector(
                                 Grid{global_ie, global_je, global_ke}) -
                                 global_grid_space->getGridOriginVector(
                                     Grid{global_is, global_js, global_ks}));
}

auto NFFFT::initNode() -> void {
  const auto grid_space = gridSpacePtr();
  setNodeBox(grid_space->getGridBoxWithoutCheck(_cube.get()));

  setNodeTaskSurfaceXN(makeNodeAxisTask(Axis::Direction::XN));
  setNodeTaskSurfaceXP(makeNodeAxisTask(Axis::Direction::XP));
  setNodeTaskSurfaceYN(makeNodeAxisTask(Axis::Direction::YN));
  setNodeTaskSurfaceYP(makeNodeAxisTask(Axis::Direction::YP));
  setNodeTaskSurfaceZN(makeNodeAxisTask(Axis::Direction::ZN));
  setNodeTaskSurfaceZP(makeNodeAxisTask(Axis::Direction::ZP));
}

auto NFFFT::makeNodeAxisTask(const Axis::Direction& direction)
    -> Divider::IndexTask {
  const auto global_box = globalBox();
  const auto node_box = nodeBox();
  const auto node_global_offset = gridSpacePtr()->globalBox().origin();
  const auto node_lower = gridSpacePtr()->box().origin();
  const auto node_upper = gridSpacePtr()->box().end();

  const auto xyz = Axis::formDirectionToXYZ(direction);

  const auto [global_offset_a, global_offset_b, global_offset_c] =
      transform::xYZToABC(
          std::tuple{node_global_offset.i(), node_global_offset.j(),
                     node_global_offset.k()},
          xyz);

  const auto [global_range_a, global_range_b, global_range_c] =
      transform::xYZToABC(
          std::tuple{Divider::makeIndexRange(global_box.origin().i(),
                                             global_box.end().i()),
                     Divider::makeIndexRange(global_box.origin().j(),
                                             global_box.end().j()),
                     Divider::makeIndexRange(global_box.origin().k(),
                                             global_box.end().k())},
          xyz);

  const auto [node_lower_a, node_lower_b, node_lower_c] = transform::xYZToABC(
      std::tuple{node_lower.i(), node_lower.j(), node_lower.k()}, xyz);

  const auto [node_upper_a, node_upper_b, node_upper_c] = transform::xYZToABC(
      std::tuple{node_upper.i(), node_upper.j(), node_upper.k()}, xyz);

  auto [node_range_a, node_range_b, node_range_c] = transform::xYZToABC(
      std::tuple{
          Divider::makeIndexRange(node_box.origin().i(), node_box.end().i()),
          Divider::makeIndexRange(node_box.origin().j(), node_box.end().j()),
          Divider::makeIndexRange(node_box.origin().k(), node_box.end().k())},
      xyz);

  auto invalid_task = Divider::IndexTask{};
  auto range_c = Divider::IndexRange{};

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

    range_c =
        Divider::makeIndexRange(node_range_c.start(), node_range_c.start() + 1);
  }

  if (node_range_c.end() + global_offset_c == global_dest_c) {
    if (node_range_c.end() >= node_upper_c) {
      // invalid
      return invalid_task;
    }

    range_c =
        Divider::makeIndexRange(node_range_c.end(), node_range_c.end() + 1);
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

    return Divider::makeIndexRange(start, end);
  };

  auto range_a = make_secondary_axis_range(
      node_lower_a, node_upper_a, node_range_a.start(), node_range_a.end());
  auto range_b = make_secondary_axis_range(
      node_lower_b, node_upper_b, node_range_b.start(), node_range_b.end());

  auto [range_x, range_y, range_z] =
      transform::aBCToXYZ(std::tuple{range_a, range_b, range_c}, xyz);

  return Divider::makeIndexTask(range_x, range_y, range_z);
}

auto NFFFT::generateSurface() -> void {
  for (const auto& f : _frequencies) {
    _fd_plane_data.emplace_back(_grid_space, _emf, f, nodeTaskSurfaceXN(),
                                nodeTaskSurfaceXP(), nodeTaskSurfaceYN(),
                                nodeTaskSurfaceYP(), nodeTaskSurfaceZN(),
                                nodeTaskSurfaceZP());
  }
}

auto NFFFT::processFarField(const xt::xtensor<double, 1>& theta,
                            const xt::xtensor<double, 1>& phi,
                            const std::string& sub_dir,
                            const Vector& origin) const -> void {
  if (!valid()) {
    return;
  }

  xt::xtensor<std::complex<double>, 1> node_data;
  for (const auto& f : _fd_plane_data) {
    const auto freq = f.frequency();
    auto node_a_theta = f.aTheta(theta, phi, origin);
    auto node_f_phi = f.fPhi(theta, phi, origin);
    auto node_a_phi = f.aPhi(theta, phi, origin);
    auto node_f_theta = f.fTheta(theta, phi, origin);

    xt::xtensor<std::complex<double>, 1> a_theta = xt::zeros_like(node_a_theta);
    xt::xtensor<std::complex<double>, 1> f_phi = xt::zeros_like(node_f_phi);
    xt::xtensor<std::complex<double>, 1> a_phi = xt::zeros_like(node_a_phi);
    xt::xtensor<std::complex<double>, 1> f_theta = xt::zeros_like(node_f_theta);

    node_data = xt::concatenate(xt::xtuple(node_data, node_a_theta));
    node_data = xt::concatenate(xt::xtuple(node_data, node_f_phi));
    node_data = xt::concatenate(xt::xtuple(node_data, node_a_phi));
    node_data = xt::concatenate(xt::xtuple(node_data, node_f_theta));
  }

  // rename later
  auto nffft_gather_func =
      [this](const xt::xtensor<std::complex<double>, 1>& send_data,
             xt::xtensor<std::complex<double>, 1>& recv_data) {
        if (this->nffftMPIConfig().size() <= 1) {
          recv_data = send_data;
          return;
        }

        MpiSupport::instance().reduceSum(this->nffftMPIConfig(),
                                         send_data.data(), recv_data.data(),
                                         send_data.size());
      };

  auto data = xt::zeros_like(node_data);
  nffft_gather_func(node_data, data);

  if (!nffftMPIConfig().isRoot()) {
    return;
  }

  const auto offsets = theta.size() * phi.size();
  for (auto i{0}; i < _fd_plane_data.size(); ++i) {
    auto freq = _fd_plane_data[i].frequency();
    std::stringstream ss;
    ss << std::fixed << std::setprecision(2);
    ss << freq / 1e9 << "GHz";
    const auto prefix_file_name = ss.str();
    if (!std::filesystem::exists(std::filesystem::path{_output_dir} /
                                 sub_dir)) {
      std::filesystem::create_directories(std::filesystem::path{_output_dir} /
                                          sub_dir);
    }

    auto&& fd_data =
        xt::view(data, xt::range(4 * i * offsets, (4 * i + 4) * offsets));
    auto a_theta = xt::view(fd_data, xt::range(0, offsets));
    auto f_phi = xt::view(fd_data, xt::range(offsets, 2 * offsets));
    auto a_phi = xt::view(fd_data, xt::range(2 * offsets, 3 * offsets));
    auto f_theta = xt::view(fd_data, xt::range(3 * offsets, 4 * offsets));

    xt::dump_npy((std::filesystem::path{_output_dir} / sub_dir /
                  (prefix_file_name + "_a_theta.npy"))
                     .string(),
                 a_theta);
    xt::dump_npy((std::filesystem::path{_output_dir} / sub_dir /
                  (prefix_file_name + "_f_phi.npy"))
                     .string(),
                 f_phi);
    xt::dump_npy((std::filesystem::path{_output_dir} / sub_dir /
                  (prefix_file_name + "_a_phi.npy"))
                     .string(),
                 a_phi);
    xt::dump_npy((std::filesystem::path{_output_dir} / sub_dir /
                  (prefix_file_name + "_f_theta.npy"))
                     .string(),
                 f_theta);
  }
}

auto NFFFT::distanceX() const -> std::size_t { return _distance_x; }

auto NFFFT::distanceY() const -> std::size_t { return _distance_y; }

auto NFFFT::distanceZ() const -> std::size_t { return _distance_z; }

auto NFFFT::cube() const -> const Cube& { return *_cube; }

auto NFFFT::globalBox() const -> const GridBox& { return _global_box; }

auto NFFFT::nodeBox() const -> const GridBox& { return _node_box; }

auto NFFFT::globalTaskSurfaceXN() const -> const Divider::IndexTask& {
  return _global_task_surface_xn;
}

auto NFFFT::globalTaskSurfaceXP() const -> const Divider::IndexTask& {
  return _global_task_surface_xp;
}

auto NFFFT::globalTaskSurfaceYN() const -> const Divider::IndexTask& {
  return _global_task_surface_yn;
}

auto NFFFT::globalTaskSurfaceYP() const -> const Divider::IndexTask& {
  return _global_task_surface_yp;
}

auto NFFFT::globalTaskSurfaceZN() const -> const Divider::IndexTask& {
  return _global_task_surface_zn;
}

auto NFFFT::globalTaskSurfaceZP() const -> const Divider::IndexTask& {
  return _global_task_surface_zp;
}

auto NFFFT::nodeTaskSurfaceXN() const -> const Divider::IndexTask& {
  return _node_task_surface_xn;
}

auto NFFFT::nodeTaskSurfaceXP() const -> const Divider::IndexTask& {
  return _node_task_surface_xp;
}

auto NFFFT::nodeTaskSurfaceYN() const -> const Divider::IndexTask& {
  return _node_task_surface_yn;
}

auto NFFFT::nodeTaskSurfaceYP() const -> const Divider::IndexTask& {
  return _node_task_surface_yp;
}

auto NFFFT::nodeTaskSurfaceZN() const -> const Divider::IndexTask& {
  return _node_task_surface_zn;
}

auto NFFFT::nodeTaskSurfaceZP() const -> const Divider::IndexTask& {
  return _node_task_surface_zp;
}

auto NFFFT::setOutputDir(const std::string& out_dir) -> void {
  _output_dir = out_dir;
}

auto NFFFT::nffftMPIConfig() -> MpiConfig& { return _nffft_mpi_config; }

auto NFFFT::nffftMPIConfig() const -> const MpiConfig& {
  return _nffft_mpi_config;
}

const GridSpace* NFFFT::gridSpacePtr() const { return _grid_space.get(); }

const CalculationParam* NFFFT::calculationParamPtr() const {
  return _calculation_param.get();
}

const EMF* NFFFT::emfPtr() const { return _emf.get(); }

auto NFFFT::setGlobalBox(GridBox box) -> void { _global_box = box; }

auto NFFFT::setNodeBox(GridBox box) -> void { _node_box = box; }

auto NFFFT::setGlobalTaskSurfaceXN(Divider::IndexTask task) -> void {
  _global_task_surface_xn = task;
}

auto NFFFT::setGlobalTaskSurfaceXP(Divider::IndexTask task) -> void {
  _global_task_surface_xp = task;
}

auto NFFFT::setGlobalTaskSurfaceYN(Divider::IndexTask task) -> void {
  _global_task_surface_yn = task;
}

auto NFFFT::setGlobalTaskSurfaceYP(Divider::IndexTask task) -> void {
  _global_task_surface_yp = task;
}

auto NFFFT::setGlobalTaskSurfaceZN(Divider::IndexTask task) -> void {
  _global_task_surface_zn = task;
}

auto NFFFT::setGlobalTaskSurfaceZP(Divider::IndexTask task) -> void {
  _global_task_surface_zp = task;
}

auto NFFFT::setNodeTaskSurfaceXN(Divider::IndexTask task) -> void {
  _node_task_surface_xn = task;
}

auto NFFFT::setNodeTaskSurfaceXP(Divider::IndexTask task) -> void {
  _node_task_surface_xp = task;
}

auto NFFFT::setNodeTaskSurfaceYN(Divider::IndexTask task) -> void {
  _node_task_surface_yn = task;
}

auto NFFFT::setNodeTaskSurfaceYP(Divider::IndexTask task) -> void {
  _node_task_surface_yp = task;
}

auto NFFFT::setNodeTaskSurfaceZN(Divider::IndexTask task) -> void {
  _node_task_surface_zn = task;
}

auto NFFFT::setNodeTaskSurfaceZP(Divider::IndexTask task) -> void {
  _node_task_surface_zp = task;
}

/* void NFFFT::initDFT() {
  using namespace std::complex_literals;
  const auto total_time_step{
      calculationParamPtr()->timeParam()->endTimeStep() -
      calculationParamPtr()->timeParam()->startTimeStep()};
  const auto dt{calculationParamPtr()->timeParam()->dt()};
  _transform_e =
      xt::zeros<std::complex<double>>({_frequencies.size(),
      total_time_step});
  _transform_h =
      xt::zeros<std::complex<double>>({_frequencies.size(),
      total_time_step});
  for (auto f{0}; f < _frequencies.size(); ++f) {
    for (std::size_t t{0}; t < total_time_step; ++t) {
      _transform_e(f, t) = dt * std::exp(-1i * 2.0 * constant::PI *
                                         (_frequencies(f) * (t + 1) * dt));
      _transform_h(f, t) = dt * std::exp(-1i * 2.0 * constant::PI *
                                         (_frequencies(f) * (t + 0.5) * dt));
    }
  }
} */

/* void NFFFT::update() {
  auto emf{emfPtr()};
  auto box{nodeBox()};
  auto is{box.origin().i()};
  auto js{box.origin().j()};
  auto ks{box.origin().k()};
  auto ie{box.end().i()};
  auto je{box.end().j()};
  auto ke{box.end().k()};
  auto current_time_step{calculationParamPtr()->timeParam()->currentTimeStep()};

  for (auto j{js}; j < je; ++j) {
    for (auto k{ks}; k < ke; ++k) {
      auto my_xn{-0.5 * (emf->ez()(is, j + 1, k) + emf->ez()(is, j, k))};
      auto mz_xn{0.5 * (emf->ey()(is, j, k + 1) + emf->ey()(is, j, k))};
      auto jy_xn{0.25 *
                 (emf->hz()(is, j, k + 1) + emf->hz()(is, j, k) +
                  emf->hz()(is - 1, j, k + 1) + emf->hz()(is - 1, j, k))};
      auto jz_xn{-0.25 *
                 (emf->hy()(is, j + 1, k) + emf->hy()(is, j, k) +
                  emf->hy()(is - 1, j + 1, k) + emf->hy()(is - 1, j, k))};

      auto my_xp{0.5 * (emf->ez()(ie, j + 1, k) + emf->ez()(ie, j, k))};
      auto mz_xp{-0.5 * (emf->ey()(ie, j, k + 1) + emf->ey()(ie, j, k))};
      auto jy_xp{-0.25 *
                 (emf->hz()(ie, j, k + 1) + emf->hz()(ie, j, k) +
                  emf->hz()(ie - 1, j, k + 1) + emf->hz()(ie - 1, j, k))};
      auto jz_xp{0.25 *
                 (emf->hy()(ie, j + 1, k) + emf->hy()(ie, j, k) +
                  emf->hy()(ie - 1, j + 1, k) + emf->hy()(ie - 1, j, k))};

      for (std::size_t f{0}; f < _frequencies.size(); ++f) {
        _my_xn(f, 0, j - js, k - ks) +=
            my_xn * _transform_e(f, current_time_step);
        _mz_xn(f, 0, j - js, k - ks) +=
            mz_xn * _transform_e(f, current_time_step);
        _jy_xn(f, 0, j - js, k - ks) +=
            jy_xn * _transform_h(f, current_time_step);
        _jz_xn(f, 0, j - js, k - ks) +=
            jz_xn * _transform_h(f, current_time_step);
        _my_xp(f, 0, j - js, k - ks) +=
            my_xp * _transform_e(f, current_time_step);
        _mz_xp(f, 0, j - js, k - ks) +=
            mz_xp * _transform_e(f, current_time_step);
        _jy_xp(f, 0, j - js, k - ks) +=
            jy_xp * _transform_h(f, current_time_step);
        _jz_xp(f, 0, j - js, k - ks) +=
            jz_xp * _transform_h(f, current_time_step);
      }
    }
  }

  for (auto i{is}; i < ie; ++i) {
    for (auto k{ks}; k < ke; ++k) {
      auto mz_yn{-0.5 * (emf->ex()(i, js, k + 1) + emf->ex()(i, js, k))};
      auto mx_yn{0.5 * (emf->ez()(i + 1, js, k) + emf->ez()(i, js, k))};
      auto jz_yn{0.25 *
                 (emf->hx()(i + 1, js, k) + emf->hx()(i, js, k) +
                  emf->hx()(i + 1, js - 1, k) + emf->hx()(i, js - 1, k))};
      auto jx_yn{-0.25 *
                 (emf->hz()(i, js, k + 1) + emf->hz()(i, js, k) +
                  emf->hz()(i, js - 1, k + 1) + emf->hz()(i, js - 1, k))};

      auto mz_yp{0.5 * (emf->ex()(i, je, k + 1) + emf->ex()(i, je, k))};
      auto mx_yp{-0.5 * (emf->ez()(i + 1, je, k) + emf->ez()(i, je, k))};
      auto jz_yp{-0.25 *
                 (emf->hx()(i + 1, je, k) + emf->hx()(i, je, k) +
                  emf->hx()(i + 1, je - 1, k) + emf->hx()(i, je - 1, k))};
      auto jx_yp{0.25 *
                 (emf->hz()(i, je, k + 1) + emf->hz()(i, je, k) +
                  emf->hz()(i, je - 1, k + 1) + emf->hz()(i, je - 1, k))};

      for (std::size_t f{0}; f < _frequencies.size(); ++f) {
        _mx_yn(f, i - is, 0, k - ks) +=
            mx_yn * _transform_e(f, current_time_step);
        _mz_yn(f, i - is, 0, k - ks) +=
            mz_yn * _transform_e(f, current_time_step);
        _jx_yn(f, i - is, 0, k - ks) +=
            jx_yn * _transform_h(f, current_time_step);
        _jz_yn(f, i - is, 0, k - ks) +=
            jz_yn * _transform_h(f, current_time_step);
        _mx_yp(f, i - is, 0, k - ks) +=
            mx_yp * _transform_e(f, current_time_step);
        _mz_yp(f, i - is, 0, k - ks) +=
            mz_yp * _transform_e(f, current_time_step);
        _jx_yp(f, i - is, 0, k - ks) +=
            jx_yp * _transform_h(f, current_time_step);
        _jz_yp(f, i - is, 0, k - ks) +=
            jz_yp * _transform_h(f, current_time_step);
      }
    }
  }

  for (auto i{is}; i < ie; ++i) {
    for (auto j{js}; j < je; ++j) {
      auto mx_zn{-0.5 * (emf->ey()(i + 1, j, ks) + emf->ey()(i, j, ks))};
      auto my_zn{0.5 * (emf->ex()(i, j + 1, ks) + emf->ex()(i, j, ks))};
      auto jx_zn{0.25 *
                 (emf->hy()(i, j + 1, ks) + emf->hy()(i, j, ks) +
                  emf->hy()(i, j + 1, ks - 1) + emf->hy()(i, j, ks - 1))};
      auto jy_zn{-0.25 *
                 (emf->hx()(i + 1, j, ks) + emf->hx()(i, j, ks) +
                  emf->hx()(i + 1, j, ks - 1) + emf->hx()(i, j, ks - 1))};

      auto mx_zp{0.5 * (emf->ey()(i + 1, j, ke) + emf->ey()(i, j, ke))};
      auto my_zp{-0.5 * (emf->ex()(i, j + 1, ke) + emf->ex()(i, j, ke))};
      auto jx_zp{-0.25 *
                 (emf->hy()(i, j + 1, ke) + emf->hy()(i, j, ke) +
                  emf->hy()(i, j + 1, ke - 1) + emf->hy()(i, j, ke - 1))};
      auto jy_zp{0.25 *
                 (emf->hx()(i + 1, j, ke) + emf->hx()(i, j, ke) +
                  emf->hx()(i + 1, j, ke - 1) + emf->hx()(i, j, ke - 1))};

      for (std::size_t f{0}; f < _frequencies.size(); ++f) {
        _mx_zn(f, i - is, j - js, 0) +=
            mx_zn * _transform_e(f, current_time_step);
        _my_zn(f, i - is, j - js, 0) +=
            my_zn * _transform_e(f, current_time_step);
        _jx_zn(f, i - is, j - js, 0) +=
            jx_zn * _transform_h(f, current_time_step);
        _jy_zn(f, i - is, j - js, 0) +=
            jy_zn * _transform_h(f, current_time_step);
        _mx_zp(f, i - is, j - js, 0) +=
            mx_zp * _transform_e(f, current_time_step);
        _my_zp(f, i - is, j - js, 0) +=
            my_zp * _transform_e(f, current_time_step);
        _jx_zp(f, i - is, j - js, 0) +=
            jx_zp * _transform_h(f, current_time_step);
        _jy_zp(f, i - is, j - js, 0) +=
            jy_zp * _transform_h(f, current_time_step);
      }
    }
  }
} */

/* void NFFFT::processFarField(const xt::xtensor<double, 1>& theta,
                            const xt::xtensor<double, 1>& phi,
                            const std::string& sub_dir, const Vector& origin) {
  using namespace std::complex_literals;
  auto box{nodeBox()};
  const auto is{box.origin().i()};
  const auto js{box.origin().j()};
  const auto ks{box.origin().k()};
  const auto ie{box.end().i()};
  const auto je{box.end().j()};
  const auto ke{box.end().k()};

  auto g{gridSpacePtr()};
  const auto& size_x{g->eSizeX()};
  const auto& size_y{g->eSizeY()};
  const auto& size_z{g->eSizeZ()};

  const auto& h_node_x{g->hNodeX()};
  const auto& h_node_y{g->hNodeY()};
  const auto& h_node_z{g->hNodeZ()};
  const auto& e_node_x{g->eNodeX()};
  const auto& e_node_y{g->eNodeY()};
  const auto& e_node_z{g->eNodeZ()};

  if (theta.size() != 1 && phi.size() != 1) {
    throw std::runtime_error("theta.size() != 1 || phi.size() != 0");
  }

  auto num_angle{theta.size() * phi.size()};

  _a_theta = xt::zeros<std::complex<double>>({_frequencies.size(), num_angle});
  _a_phi = xt::zeros<std::complex<double>>({_frequencies.size(), num_angle});
  _f_theta = xt::zeros<std::complex<double>>({_frequencies.size(), num_angle});
  _f_phi = xt::zeros<std::complex<double>>({_frequencies.size(), num_angle});

  auto cos_t{xt::cos(theta)};
  auto sin_t{xt::sin(theta)};
  auto cos_p{xt::cos(phi)};
  auto sin_p{xt::sin(phi)};
  auto sin_t_sin_p{sin_t * sin_p};
  auto sin_t_cos_p{sin_t * cos_p};
  auto cos_t_sin_p{cos_t * sin_p};
  auto cos_t_cos_p{cos_t * cos_p};

  for (std::size_t f{0}; f < _frequencies.size(); ++f) {
    auto wave_number{2.0 * constant::PI * _frequencies(f) / constant::C_0};

    for (auto j{js}; j < je; ++j) {
      for (auto k{ks}; k < ke; ++k) {
        auto ds{size_y(j) * size_z(k)};
        auto r_xn{Vector{e_node_x(is), h_node_y(j), h_node_z(k)} - origin};
        auto phase_shift_xn{
            xt::exp(1i * wave_number *
                    (r_xn.x() * sin_t_cos_p + r_xn.y() * sin_t_sin_p +
                     r_xn.z() * cos_t))};
        xt::view(_a_theta, f, xt::all()) +=
            (_jy_xn(f, 0, j - js, k - ks) * cos_t_sin_p -
             _jz_xn(f, 0, j - js, k - ks) * sin_t) *
            phase_shift_xn * ds;
        xt::view(_f_theta, f, xt::all()) +=
            (_my_xn(f, 0, j - js, k - ks) * cos_t_sin_p -
             _mz_xn(f, 0, j - js, k - ks) * sin_t) *
            phase_shift_xn * ds;
        xt::view(_a_phi, f, xt::all()) +=
            (_jy_xn(f, 0, j - js, k - ks) * cos_p) * phase_shift_xn * ds;
        xt::view(_f_phi, f, xt::all()) +=
            (_my_xn(f, 0, j - js, k - ks) * cos_p) * phase_shift_xn * ds;
        auto r_xp{Vector{e_node_x(ie), h_node_y(j), h_node_z(k)} - origin};
        auto phase_shift_xp{
            xt::exp(1i * wave_number *
                    (r_xp.x() * sin_t_cos_p + r_xp.y() * sin_t_sin_p +
                     r_xp.z() * cos_t))};
        xt::view(_a_theta, f, xt::all()) +=
            (_jy_xp(f, 0, j - js, k - ks) * cos_t_sin_p -
             _jz_xp(f, 0, j - js, k - ks) * sin_t) *
            phase_shift_xp * ds;
        xt::view(_f_theta, f, xt::all()) +=
            (_my_xp(f, 0, j - js, k - ks) * cos_t_sin_p -
             _mz_xp(f, 0, j - js, k - ks) * sin_t) *
            phase_shift_xp * ds;
        xt::view(_a_phi, f, xt::all()) +=
            (_jy_xp(f, 0, j - js, k - ks) * cos_p) * phase_shift_xp * ds;
        xt::view(_f_phi, f, xt::all()) +=
            (_my_xp(f, 0, j - js, k - ks) * cos_p) * phase_shift_xp * ds;
      }
    }

    for (auto i{is}; i < ie; ++i) {
      for (auto k{ks}; k < ke; ++k) {
        auto ds{size_x(i) * size_z(k)};
        auto r_yn{Vector{h_node_x(i), e_node_y(js), h_node_z(k)} - origin};
        auto phase_shift_yn{
            xt::exp(1i * wave_number *
                    (r_yn.x() * sin_t_cos_p + r_yn.y() * sin_t_sin_p +
                     r_yn.z() * cos_t))};
        xt::view(_a_theta, f, xt::all()) +=
            (_jx_yn(f, i - is, 0, k - ks) * cos_t_cos_p -
             _jz_yn(f, i - is, 0, k - ks) * sin_t) *
            phase_shift_yn * ds;
        xt::view(_f_theta, f, xt::all()) +=
            (_mx_yn(f, i - is, 0, k - ks) * cos_t_cos_p -
             _mz_yn(f, i - is, 0, k - ks) * sin_t) *
            phase_shift_yn * ds;
        xt::view(_a_phi, f, xt::all()) -=
            (_jx_yn(f, i - is, 0, k - ks) * sin_p) * phase_shift_yn * ds;
        xt::view(_f_phi, f, xt::all()) -=
            (_mx_yn(f, i - is, 0, k - ks) * sin_p) * phase_shift_yn * ds;
        auto r_yp{Vector{h_node_x(i), e_node_y(je), h_node_z(k)} - origin};
        auto phase_shift_yp{
            xt::exp(1i * wave_number *
                    (r_yp.x() * sin_t_cos_p + r_yp.y() * sin_t_sin_p +
                     r_yp.z() * cos_t))};
        xt::view(_a_theta, f, xt::all()) +=
            (_jx_yp(f, i - is, 0, k - ks) * cos_t_cos_p -
             _jz_yp(f, i - is, 0, k - ks) * sin_t) *
            phase_shift_yp * ds;
        xt::view(_f_theta, f, xt::all()) +=
            (_mx_yp(f, i - is, 0, k - ks) * cos_t_cos_p -
             _mz_yp(f, i - is, 0, k - ks) * sin_t) *
            phase_shift_yp * ds;
        xt::view(_a_phi, f, xt::all()) -=
            (_jx_yp(f, i - is, 0, k - ks) * sin_p) * phase_shift_yp * ds;
        xt::view(_f_phi, f, xt::all()) -=
            (_mx_yp(f, i - is, 0, k - ks) * sin_p) * phase_shift_yp * ds;
      }
    }

    for (auto i{is}; i < ie; ++i) {
      for (auto j{js}; j < je; ++j) {
        auto ds = size_x(i) * size_y(j);
        auto r_zn{Vector{h_node_x(i), h_node_y(j), e_node_z(ks)} - origin};
        auto phase_shift_zn{
            xt::exp(1i * wave_number *
                    (r_zn.x() * sin_t_cos_p + r_zn.y() * sin_t_sin_p +
                     r_zn.z() * cos_t))};
        xt::view(_a_theta, f, xt::all()) +=
            (_jx_zn(f, i - is, j - js, 0) * cos_t_cos_p +
             _jy_zn(f, i - is, j - js, 0) * cos_t_sin_p) *
            phase_shift_zn * ds;
        xt::view(_f_theta, f, xt::all()) +=
            (_mx_zn(f, i - is, j - js, 0) * cos_t_cos_p +
             _my_zn(f, i - is, j - js, 0) * cos_t_sin_p) *
            phase_shift_zn * ds;
        xt::view(_a_phi, f, xt::all()) +=
            (-_jx_zn(f, i - is, j - js, 0) * sin_p +
             _jy_zn(f, i - is, j - js, 0) * cos_p) *
            phase_shift_zn * ds;
        xt::view(_f_phi, f, xt::all()) +=
            (-_mx_zn(f, i - is, j - js, 0) * sin_p +
             _my_zn(f, i - is, j - js, 0) * cos_p) *
            phase_shift_zn * ds;
        auto r_zp{Vector{h_node_x(i), h_node_y(j), e_node_z(ke)} - origin};
        auto phase_shift_zp{
            xt::exp(1i * wave_number *
                    (r_zp.x() * sin_t_cos_p + r_zp.y() * sin_t_sin_p +
                     r_zp.z() * cos_t))};
        xt::view(_a_theta, f, xt::all()) +=
            (_jx_zp(f, i - is, j - js, 0) * cos_t_cos_p +
             _jy_zp(f, i - is, j - js, 0) * cos_t_sin_p) *
            phase_shift_zp * ds;
        xt::view(_f_theta, f, xt::all()) +=
            (_mx_zp(f, i - is, j - js, 0) * cos_t_cos_p +
             _my_zp(f, i - is, j - js, 0) * cos_t_sin_p) *
            phase_shift_zp * ds;
        xt::view(_a_phi, f, xt::all()) +=
            (-_jx_zp(f, i - is, j - js, 0) * sin_p +
             _jy_zp(f, i - is, j - js, 0) * cos_p) *
            phase_shift_zp * ds;
        xt::view(_f_phi, f, xt::all()) +=
            (-_mx_zp(f, i - is, j - js, 0) * sin_p +
             _my_zp(f, i - is, j - js, 0) * cos_p) *
            phase_shift_zp * ds;
      }
    }
  }

  auto out_dir{std::filesystem::path{_output_dir} / sub_dir};
  if (!std::filesystem::exists(out_dir) ||
      !std::filesystem::is_directory(out_dir)) {
    std::filesystem::create_directories(out_dir);
  }

  auto a_theta_path{out_dir / "a_theta.npy"};
  auto a_phi_path{out_dir / "a_phi.npy"};
  auto f_theta_path{out_dir / "f_theta.npy"};
  auto f_phi_path{out_dir / "f_phi.npy"};

  xt::dump_npy(a_theta_path.string(), _a_theta);
  xt::dump_npy(a_phi_path.string(), _a_phi);
  xt::dump_npy(f_theta_path.string(), _f_theta);
  xt::dump_npy(f_phi_path.string(), _f_phi);
} */

/* void NFFFT::outputRadiationPower() {
  if (!valid()) {
    return;
  }

  using namespace std::complex_literals;

  auto g{gridSpacePtr()};
  const auto& size_x{_grid_space->eSizeX()};
  const auto& size_y{_grid_space->eSizeY()};
  const auto& size_z{_grid_space->eSizeZ()};

  auto calculate_power =
      [](const auto& num_freq, const auto& direction, const auto& task,
         const auto& size_x, const auto& size_y, const auto& size_z,
         const auto& ja, const auto& jb, const auto& ma, const auto& mb) {
        xt::xarray<std::complex<double>> power =
            xt::zeros<std::complex<double>>({num_freq});

        if (!task.valid()) {
          return power;
        }
        const auto is = task.xRange().start();
        const auto ie = task.xRange().end();
        const auto js = task.yRange().start();
        const auto je = task.yRange().end();
        const auto ks = task.zRange().start();
        const auto ke = task.zRange().end();
        for (std::size_t f{0}; f < num_freq; ++f) {
          for (auto i{is}; i < ie; ++i) {
            for (auto j{js}; j < je; ++j) {
              for (auto k{ks}; k < ke; ++k) {
                power(f) += size_x(i) * size_y(j) * size_z(k) *
                            (mb.at(f, i - is, j - js, k - ks) *
                                 std::conj(ja.at(f, i - is, j - js, k - ks))
                                 -
                             ma.at(f, i - is, j - js, k - ks) *
                                 std::conj(jb.at(f, i - is, j - js, k -
                                 ks)));
              }
            }
          }
        }
        if (Axis::directionNegative(direction)) {
          power *= -1.0;
        }
        return power;
      };
  struct AlwaysOne {
    auto operator()(std::size_t i) const -> auto { return 1; }
  } always_one;

  auto future_arr =
      std::vector<std::future<xt::xarray<std::complex<double>>>>();

  future_arr.emplace_back(std::async(calculate_power, _frequencies.size(),
                                     Axis::Direction::XN,
                                     nodeTaskSurfaceXN(), always_one, size_y,
                                     size_z, _jy_xn, _jz_xn, _my_xn,
                                     _mz_xn));
  future_arr.emplace_back(std::async(calculate_power, _frequencies.size(),
                                     Axis::Direction::XP,
                                     nodeTaskSurfaceXP(), always_one, size_y,
                                     size_z, _jy_xp, _jz_xp, _my_xp,
                                     _mz_xp));
  future_arr.emplace_back(std::async(calculate_power, _frequencies.size(),
                                     Axis::Direction::YN,
                                     nodeTaskSurfaceYN(), size_x, always_one,
                                     size_z, _jz_yn, _jx_yn, _mz_yn,
                                     _mx_yn));
  future_arr.emplace_back(std::async(calculate_power, _frequencies.size(),
                                     Axis::Direction::YP,
                                     nodeTaskSurfaceYP(), size_x, always_one,
                                     size_z, _jz_yp, _jx_yp, _mz_yp,
                                     _mx_yp));
  future_arr.emplace_back(
      std::async(std::launch::deferred, calculate_power, _frequencies.size(),
                 Axis::Direction::ZN, nodeTaskSurfaceZN(), size_x, size_y,
                 always_one, _jx_zn, _jy_zn, _mx_zn, _my_zn));
  future_arr.emplace_back(
      std::async(std::launch::deferred, calculate_power, _frequencies.size(),
                 Axis::Direction::ZP, nodeTaskSurfaceZP(), size_x, size_y,
                 always_one, _jx_zp, _jy_zp, _mx_zp, _my_zp));

  xt::xarray<std::complex<double>> power =
      xt::zeros<std::complex<double>>({_frequencies.size()});

  std::for_each(future_arr.begin(), future_arr.end(), [&power](auto& future)
  {
    auto res = future.get();
    // std::cout << "RES: " << res << "\n";
    power += res;
  });

  xt::xarray<double> node_res = 0.5 * xt::real(power);
  auto global_res = xt::empty_like(node_res);

  auto gather_power = [this](const auto& node_data, auto&& data) {
    if (!this->valid()) {
      return;
    }

    auto& mpi_support = MpiSupport::instance();
    const auto& mpi_config = this->nffftMPIConfig();

    if (mpi_config.size() <= 1) {
      data = node_data;
      return;
    }

    mpi_support.reduceSum(this->nffftMPIConfig(), node_data.data(),
    data.data(),
                          data.size());
  };

  auto collect_func = [this](const auto& tag, auto&& gather_func,
                             const auto& node_data, auto&& data) {
    auto& mpi_support = MpiSupport::instance();
    const auto& mpi_config = this->nffftMPIConfig();
    const auto& node_res = xt::xarray<double>{};

    if (!this->valid()) {
      if (!mpi_support.isRoot()) {
        return;
      }

      if (mpi_support.size() <= 1) {
        throw XFDTDNFFFTException(
            "XFDTD NFFFT Exception: NFFFT is not valid, but MPI size is not "
            "1.");
      }

      mpi_support.recv(mpi_support.config(), data.data(),
                       sizeof(double) * data.size(), MpiSupport::ANY_SOURCE,
                       tag);
      return;
    }

    gather_func(node_data, data);

    if (mpi_support.size() <= 1) {
      return;
    }

    if (!mpi_support.isRoot() && mpi_config.isRoot()) {
      mpi_support.send(mpi_support.config(), data.data(),
                       sizeof(double) * data.size(),
                       mpi_support.config().root(), tag);
      return;
    }

    if (mpi_support.isRoot() && !mpi_config.isRoot()) {
      mpi_support.recv(mpi_support.config(), data.data(),
                       sizeof(double) * data.size(), MpiSupport::ANY_SOURCE,
                       tag);
      return;
    }
  };

  collect_func(41, gather_power, node_res, global_res);

  if (!MpiSupport::instance().isRoot()) {
    return;
  }

  auto radiation_power_path{std::filesystem::path{_output_dir} /
                            "radiation_power.npy"};
  if (!std::filesystem::exists(radiation_power_path) ||
      !std::filesystem::is_regular_file(radiation_power_path)) {
    std::filesystem::create_directories(radiation_power_path.parent_path());
  }
  xt::dump_npy(radiation_power_path.string(), global_res);
} */

}  // namespace xfdtd
