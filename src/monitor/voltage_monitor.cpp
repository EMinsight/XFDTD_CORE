#include <xfdtd/coordinate_system/coordinate_system.h>
#include <xfdtd/grid_space/grid_space.h>
#include <xfdtd/monitor/monitor.h>
#include <xfdtd/monitor/voltage_monitor.h>
#include <xfdtd/parallel/mpi_support.h>
#include <xfdtd/util/transform.h>

#include <sstream>
#include <tuple>
#include <utility>
#include <xtensor.hpp>

namespace xfdtd {

VoltageMonitor::VoltageMonitor(std::string name, std::unique_ptr<Cube> cube,
                               Axis::Direction direction,
                               std::string output_dir)
    : TimeMonitor{std::move(name), std::move(cube), std::move(output_dir)},
      _direction{direction} {}

void VoltageMonitor::init(
    std::shared_ptr<const GridSpace> grid_space,
    std::shared_ptr<const CalculationParam> calculation_param,
    std::shared_ptr<const EMF> emf) {
  defaultInit(std::move(grid_space), std::move(calculation_param),
              std::move(emf));

  auto new_range = [](const auto& node_lower, const auto& node_upper,
                      const IndexRange& range) {
    auto s = range.start();
    auto e = range.end();
    if (s <= node_lower) {
      s = node_lower + 1;
    }
    if (node_upper <= e) {
      e = node_upper - 1;
    }
    return makeIndexRange(s, e);
  };

  auto correct_global_abc_range =
      [&new_range](const auto& a_lower, const auto& a_upper,
                   IndexRange& a_range, const auto& b_lower,
                   const auto& b_upper, IndexRange& b_range) {
        a_range = new_range(a_lower, a_upper, a_range);
        b_range = new_range(b_lower, b_upper, b_range);
      };

  auto correct_node_abc_range =
      [&new_range](const auto& a_lower, const auto& a_upper,
                   IndexRange& a_range, const auto& b_lower,
                   const auto& b_upper, IndexRange& b_range,
                   const auto& global_c_lower, const auto global_c_upper,
                   const auto& offset_c, const auto& c_lower,
                   const auto& c_upper, IndexRange& c_range) {
        a_range = new_range(a_lower, a_upper, a_range);
        b_range = new_range(b_lower, b_upper, b_range);
        auto c_s = c_range.start();
        auto c_e = c_range.end();
        c_range = new_range(c_lower, c_upper, c_range);
        if (c_s + offset_c == global_c_lower) {
          c_range = makeIndexRange(c_s, c_range.end());
        }
        if (c_e + offset_c == global_c_upper) {
          c_range = makeIndexRange(c_range.start(), c_e);
        }
      };

  const auto global_offset = gridSpacePtr()->globalBox().origin();
  const auto [global_offset_a, global_offset_b, global_offset_c] =
      xfdtd::transform::xYZToABC(
          std::tuple(global_offset.i(), global_offset.j(), global_offset.k()),
          Axis::fromDirectionToXYZ(_direction));

  const auto global_lower = gridSpacePtr()->globalGridSpace()->box().origin();
  const auto [global_lower_a, global_lower_b, global_lower_c] =
      xfdtd::transform::xYZToABC(
          std::tuple(global_lower.i(), global_lower.j(), global_lower.k()),
          Axis::fromDirectionToXYZ(_direction));

  const auto global_upper = gridSpacePtr()->globalGridSpace()->box().end();
  const auto [global_upper_a, global_upper_b, global_upper_c] =
      xfdtd::transform::xYZToABC(
          std::tuple(global_upper.i(), global_upper.j(), global_upper.k()),
          Axis::fromDirectionToXYZ(_direction));

  // correct global task
  auto g_abc = xfdtd::transform::xYZToABC(
      std::tuple(globalTask().xRange(), globalTask().yRange(),
                 globalTask().zRange()),
      Axis::fromDirectionToXYZ(_direction));
  auto global_range_a =
      makeIndexRange(std::get<0>(g_abc).start(), std::get<0>(g_abc).end() + 1);
  auto global_range_b =
      makeIndexRange(std::get<1>(g_abc).start(), std::get<1>(g_abc).end() + 1);
  auto global_range_c =
      makeIndexRange(std::get<2>(g_abc).start(), std::get<2>(g_abc).end());
  correct_global_abc_range(global_lower_a, global_upper_a, global_range_a,
                           global_lower_b, global_upper_b, global_range_b);
  auto [global_range_x, global_range_y, global_range_z] =
      xfdtd::transform::aBCToXYZ(
          std::tuple(global_range_a, global_range_b, global_range_c),
          Axis::fromDirectionToXYZ(_direction));
  setGlobalTask(makeIndexTask(global_range_x, global_range_y, global_range_z));
  if (!globalTask().valid()) {
    throw XFDTDMonitorException(name() + " globalTask is not valid\n" +
                                globalTask().toString());
  }

  // correct node task
  const auto node_lower = gridSpacePtr()->box().origin();
  const auto [node_lower_a, node_lower_b, node_lower_c] =
      xfdtd::transform::xYZToABC(
          std::tuple(node_lower.i(), node_lower.j(), node_lower.k()),
          Axis::fromDirectionToXYZ(_direction));
  const auto node_upper = gridSpacePtr()->box().end();
  const auto [node_upper_a, node_upper_b, node_upper_c] =
      xfdtd::transform::xYZToABC(
          std::tuple(node_upper.i(), node_upper.j(), node_upper.k()),
          Axis::fromDirectionToXYZ(_direction));

  auto node_abc = xfdtd::transform::xYZToABC(
      std::tuple(nodeTask().xRange(), nodeTask().yRange(), nodeTask().zRange()),
      Axis::fromDirectionToXYZ(_direction));
  auto node_range_a = makeIndexRange(std::get<0>(node_abc).start(),
                                     std::get<0>(node_abc).end() + 1);
  auto node_range_b = makeIndexRange(std::get<1>(node_abc).start(),
                                     std::get<1>(node_abc).end() + 1);
  auto node_range_c = makeIndexRange(std::get<2>(node_abc).start(),
                                     std::get<2>(node_abc).end());

  correct_node_abc_range(node_lower_a, node_upper_a, node_range_a, node_lower_b,
                         node_upper_b, node_range_b, global_range_c.start(),
                         global_range_c.end(), global_offset_c, node_lower_c,
                         node_upper_c, node_range_c);

  auto [node_x, node_y, node_z] = xfdtd::transform::aBCToXYZ(
      std::tuple(node_range_a, node_range_b, node_range_c),
      Axis::fromDirectionToXYZ(_direction));

  setNodeTask(makeIndexTask(node_x, node_y, node_z));

  if (!nodeTask().valid()) {
    return;
  }

  // init for recording
  _is = nodeTask().xRange().start();
  _ie = nodeTask().xRange().end();
  _js = nodeTask().yRange().start();
  _je = nodeTask().yRange().end();
  _ks = nodeTask().zRange().start();
  _ke = nodeTask().zRange().end();

  if (Axis::fromDirectionToXYZ(_direction) == Axis::XYZ::X) {
    _dc = xt::view(gridSpacePtr()->hSizeX(), xt::range(_is, _ie));

    auto global_js = globalTask().yRange().start();
    auto global_je = globalTask().yRange().end();
    auto global_ks = globalTask().zRange().start();
    auto global_ke = globalTask().zRange().end();

    _coff = -_dc / ((global_je - global_js) * (global_ke - global_ks));
    if (_direction == Axis::Direction::XN) {
      _coff *= -1.0;
    }
  }

  if (Axis::fromDirectionToXYZ(_direction) == Axis::XYZ::Y) {
    _dc = xt::view(gridSpacePtr()->hSizeY(), xt::range(_js, _je));

    auto global_ks = globalTask().zRange().start();
    auto global_ke = globalTask().zRange().end();
    auto global_is = globalTask().xRange().start();
    auto global_ie = globalTask().xRange().end();

    _coff = -_dc / ((global_ke - global_ks) * (global_ie - global_is));
    if (_direction == Axis::Direction::YN) {
      _coff *= -1.0;
    }
  }

  if (Axis::fromDirectionToXYZ(_direction) == Axis::XYZ::Z) {
    _dc = xt::view(gridSpacePtr()->hSizeZ(), xt::range(_ks, _ke));

    auto global_is = globalTask().xRange().start();
    auto global_ie = globalTask().xRange().end();
    auto global_js = globalTask().yRange().start();
    auto global_je = globalTask().yRange().end();

    _coff = -_dc / ((global_ie - global_is) * (global_je - global_js));
    if (Axis::directionNegative(_direction)) {
      _coff *= -1.0;
    }
  }
}

void VoltageMonitor::update() {
  if (!valid()) {
    return;
  }

  auto emf{emfPtr()};
  auto t{calculationParamPtr()->timeParam()->currentTimeStep()};

  if (Axis::fromDirectionToXYZ(_direction) == Axis::XYZ::X) {
    for (auto i{_is}; i < _ie; ++i) {
      for (auto j{_js}; j < _je; ++j) {
        for (auto k{_ks}; k < _ke; ++k) {
          _node_data(t) += _coff(i - _is) * emf->ex()(i, j, k);
        }
      }
    }
  }

  if (Axis::fromDirectionToXYZ(_direction) == Axis::XYZ::Y) {
    for (auto i{_is}; i < _ie; ++i) {
      for (auto j{_js}; j < _je; ++j) {
        for (auto k{_ks}; k < _ke; ++k) {
          _node_data(t) += _coff(j - _js) * emf->ey()(i, j, k);
        }
      }
    }
  }

  if (Axis::fromDirectionToXYZ(_direction) == Axis::XYZ::Z) {
    for (auto i{_is}; i < _ie; ++i) {
      for (auto j{_js}; j < _je; ++j) {
        for (auto k{_ks}; k < _ke; ++k) {
          _node_data(t) += _coff(k - _ks) * emf->ez()(i, j, k);
        }
      }
    }
  }
}

void VoltageMonitor::initTimeDependentVariable() {
  setTime(calculationParamPtr()->timeParam()->eTime());
  _node_data = xt::zeros_like(time());
  data() = xt::zeros_like(time());
}

void VoltageMonitor::initParallelizedConfig() { makeMpiSubComm(); }

void VoltageMonitor::gatherData() {
  if (!valid()) {
    return;
  }

  if (monitorMpiConfig().size() <= 1) {
    data() = _node_data;
    return;
  }

  auto& mpi_support = MpiSupport::instance();

  data() = xt::zeros_like(_node_data);

  mpi_support.reduceSum(monitorMpiConfig(), _node_data.data(), data().data(),
                        data().size());
}

auto VoltageMonitor::toString() const -> std::string {
  std::stringstream ss;
  ss << "Voltage Monitor:\n";
  ss << " Global " << globalGridBox().toString() << "\n";
  ss << " Global " << globalTask().toString() << "\n";
  ss << " Node " << nodeTask().toString() << "\n";
  return ss.str();
}

}  // namespace xfdtd
