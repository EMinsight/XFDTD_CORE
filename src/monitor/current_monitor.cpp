#include <xfdtd/coordinate_system/coordinate_system.h>
#include <xfdtd/divider/divider.h>
#include <xfdtd/monitor/current_monitor.h>
#include <xfdtd/monitor/monitor.h>
#include <xfdtd/parallel/mpi_support.h>
#include <xfdtd/util/transform.h>

#include <memory>
#include <sstream>
#include <string>
#include <utility>
#include <xtensor.hpp>

namespace xfdtd {

CurrentMonitor::CurrentMonitor(std::string name, std::unique_ptr<Cube> cube,
                               Axis::Direction direction,
                               std::string output_dir)
    : TimeMonitor{std::move(name), std::move(cube), std::move(output_dir)},
      _direction{direction} {}

void CurrentMonitor::init(
    std::shared_ptr<const GridSpace> grid_space,
    std::shared_ptr<const CalculationParam> calculation_param,
    std::shared_ptr<const EMF> emf) {
  defaultInit(grid_space, calculation_param, emf);

  _is = nodeGridBox().origin().i();
  _ie = nodeGridBox().end().i();
  _js = nodeGridBox().origin().j();
  _je = nodeGridBox().end().j();
  _ks = nodeGridBox().origin().k();
  _ke = nodeGridBox().end().k();

  const auto global_offset = gridSpacePtr()->globalBox().origin();
  const auto [global_offset_a, global_offset_b, global_offset_c] =
      xfdtd::transform::xYZToABC(
          std::tuple(global_offset.i(), global_offset.j(), global_offset.k()),
          Axis::formDirectionToXYZ(_direction));

  const auto global_lower = gridSpacePtr()->globalGridSpace()->box().origin();
  const auto [global_lower_a, global_lower_b, global_lower_c] =
      xfdtd::transform::xYZToABC(
          std::tuple(global_lower.i(), global_lower.j(), global_lower.k()),
          Axis::formDirectionToXYZ(_direction));

  const auto global_upper = gridSpacePtr()->globalGridSpace()->box().end();
  const auto [global_upper_a, global_upper_b, global_upper_c] =
      xfdtd::transform::xYZToABC(
          std::tuple(global_upper.i(), global_upper.j(), global_upper.k()),
          Axis::formDirectionToXYZ(_direction));

  // correct global task
  auto g_abc = xfdtd::transform::xYZToABC(
      std::tuple(globalTask().xRange(), globalTask().yRange(),
                 globalTask().zRange()),
      Axis::formDirectionToXYZ(_direction));
  auto global_range_c = Divider::makeIndexRange(std::get<2>(g_abc).end() - 1,
                                                std::get<2>(g_abc).end());
  auto global_c_e = global_range_c.end();
  if (global_c_e == 0) {
    throw XFDTDMonitorException("CurrentMonitor::init(): global_c_e is 0");
  }

  if (std::get<0>(g_abc).start() == 0 ||
      std::get<0>(g_abc).end() == global_upper_a) {
    throw XFDTDMonitorException(
        std::string("CurrentMonitor::init(): global x range is invalid: ") +
        std::get<0>(g_abc).toString());
  }

  if (std::get<1>(g_abc).start() == 0 ||
      std::get<1>(g_abc).end() == global_upper_b) {
    throw XFDTDMonitorException(
        std::string("CurrentMonitor::init(): global y range is invalid: ") +
        std::get<1>(g_abc).toString());
  }

  auto [global_range_x, global_range_y, global_range_z] =
      xfdtd::transform::aBCToXYZ(
          std::tuple(std::get<0>(g_abc), std::get<1>(g_abc), global_range_c),
          Axis::formDirectionToXYZ(_direction));
  setGlobalTask(
      Divider::makeIndexTask(global_range_x, global_range_y, global_range_z));

  // correct node task
  const auto node_lower = gridSpacePtr()->box().origin();
  const auto [node_lower_a, node_lower_b, node_lower_c] =
      xfdtd::transform::xYZToABC(
          std::tuple(node_lower.i(), node_lower.j(), node_lower.k()),
          Axis::formDirectionToXYZ(_direction));
  const auto node_upper = gridSpacePtr()->box().end();
  const auto [node_upper_a, node_upper_b, node_upper_c] =
      xfdtd::transform::xYZToABC(
          std::tuple(node_upper.i(), node_upper.j(), node_upper.k()),
          Axis::formDirectionToXYZ(_direction));

  auto n_abc = xfdtd::transform::xYZToABC(
      std::tuple(nodeTask().xRange(), nodeTask().yRange(), nodeTask().zRange()),
      Axis::formDirectionToXYZ(_direction));
  auto node_range_a = Divider::makeIndexRange(std::get<0>(n_abc).start(),
                                              std::get<0>(n_abc).end());
  auto node_range_b = Divider::makeIndexRange(std::get<1>(n_abc).start(),
                                              std::get<1>(n_abc).end());
  auto node_range_c = Divider::makeIndexRange(std::get<2>(n_abc).end() - 1,
                                              std::get<2>(n_abc).end());

  if (node_range_c.end() == 0) {
    return;
  }

  if (node_range_c.end() + global_offset_c != global_c_e) {
    // make it invalid
    setNodeTask(Divider::makeIndexTask(
        nodeTask().xRange(), nodeTask().yRange(),
        Divider::makeIndexRange(node_range_c.end(), node_range_c.end())));
    return;
  }

  auto [node_range_x, node_range_y, node_range_z] = xfdtd::transform::aBCToXYZ(
      std::tuple(node_range_a, node_range_b, node_range_c),
      Axis::formDirectionToXYZ(_direction));
  setNodeTask(Divider::makeIndexTask(node_range_x, node_range_y, node_range_z));
  _is = node_range_x.start();
  _ie = node_range_x.end();
  _js = node_range_y.start();
  _je = node_range_y.end();
  _ks = node_range_z.start();
  _ke = node_range_z.end();

  auto make_ranges = [](const auto& node_task_b_start,
                        const auto& node_task_b_end, const auto& node_b_upper,
                        const auto& global_b_offset, const auto& global_b_start,
                        const auto& global_b_end, const auto& node_a_lower,
                        const auto& node_a_upper, auto start, auto end) {
    if (start <= node_a_lower) {
      start = node_a_lower + 1;
    }
    end = end + 1;
    if (node_a_upper <= end) {
      end = node_a_upper - 1;
    }

    auto range_n = Divider::makeIndexRange(0, 0);
    auto range_p = Divider::makeIndexRange(0, 0);

    if (node_task_b_start != 0 && node_task_b_start != node_b_upper - 1 &&
        node_task_b_start != node_b_upper &&
        node_task_b_start + global_b_offset == global_b_start) {
      range_n = Divider::makeIndexRange(start, end);
    }

    if (node_task_b_end != node_b_upper &&
        node_task_b_end != node_b_upper - 1 && node_task_b_end != 0 &&
        node_task_b_end + global_b_offset == global_b_end) {
      range_p = Divider::makeIndexRange(start, end);
    }

    return std::pair{range_n, range_p};
  };

  auto [ha_bn, ha_bp] = make_ranges(
      node_range_b.start(), node_range_b.end(), node_upper_b, global_offset_b,
      std::get<1>(g_abc).start(), std::get<1>(g_abc).end(), node_lower_a,
      node_upper_a, node_range_a.start(), node_range_a.end());

  _ha_range_bn = ha_bn;
  _ha_range_bp = ha_bp;

  auto [hb_an, hb_ap] = make_ranges(
      node_range_a.start(), node_range_a.end(), node_upper_a, global_offset_a,
      std::get<0>(g_abc).start(), std::get<0>(g_abc).end(), node_lower_b,
      node_upper_b, node_range_b.start(), node_range_b.end());

  _hb_range_an = hb_an;
  _hb_range_ap = hb_ap;

  if (Axis::formDirectionToXYZ(_direction) == Axis::XYZ::X) {
    // _da = xt::view(gridSpacePtr()->eSizeY(), xt::range(_js, _je + 1));
    // _db = xt::view(gridSpacePtr()->eSizeZ(), xt::range(_ks, _ke + 1));
    if (_ha_range_bn.valid()) {
      _da = xt::view(gridSpacePtr()->eSizeY(),
                     xt::range(_ha_range_bn.start(), _ha_range_bn.end()));
    } else if (_ha_range_bp.valid()) {
      _da = xt::view(gridSpacePtr()->eSizeY(),
                     xt::range(_ha_range_bp.start(), _ha_range_bp.end()));
    }

    if (_hb_range_an.valid()) {
      _db = xt::view(gridSpacePtr()->eSizeZ(),
                     xt::range(_hb_range_an.start(), _hb_range_an.end()));
    } else if (_hb_range_ap.valid()) {
      _db = xt::view(gridSpacePtr()->eSizeZ(),
                     xt::range(_hb_range_ap.start(), _hb_range_ap.end()));
    }
    if (Axis::directionPositive(_direction)) {
      _positive = 1.0;
    } else {
      _positive = -1.0;
    }
  }

  if (Axis::formDirectionToXYZ(_direction) == Axis::XYZ::Y) {
    // _da = xt::view(gridSpacePtr()->eSizeZ(), xt::range(_ks, _ke + 1));
    // _db = xt::view(gridSpacePtr()->eSizeX(), xt::range(_is, _ie + 1));
    if (_ha_range_bn.valid()) {
      _da = xt::view(gridSpacePtr()->eSizeZ(),
                     xt::range(_ha_range_bn.start(), _ha_range_bn.end()));
    } else if (_ha_range_bp.valid()) {
      _da = xt::view(gridSpacePtr()->eSizeZ(),
                     xt::range(_ha_range_bp.start(), _ha_range_bp.end()));
    }

    if (_hb_range_an.valid()) {
      _db = xt::view(gridSpacePtr()->eSizeX(),
                     xt::range(_hb_range_an.start(), _hb_range_an.end()));
    } else if (_hb_range_ap.valid()) {
      _db = xt::view(gridSpacePtr()->eSizeX(),
                     xt::range(_hb_range_ap.start(), _hb_range_ap.end()));
    }
    if (Axis::directionPositive(_direction)) {
      _positive = 1.0;
    } else {
      _positive = -1.0;
    }
  }

  if (Axis::formDirectionToXYZ(_direction) == Axis::XYZ::Z) {
    // _da = xt::view(gridSpacePtr()->eSizeX(), xt::range(_is, _ie + 1));
    // _db = xt::view(gridSpacePtr()->eSizeY(), xt::range(_js, _je + 1));
    if (_ha_range_bn.valid()) {
      _da = xt::view(gridSpacePtr()->eSizeX(),
                     xt::range(_ha_range_bn.start(), _ha_range_bn.end()));
    } else if (_ha_range_bp.valid()) {
      _da = xt::view(gridSpacePtr()->eSizeX(),
                     xt::range(_ha_range_bp.start(), _ha_range_bp.end()));
    }

    if (_hb_range_an.valid()) {
      _db = xt::view(gridSpacePtr()->eSizeY(),
                     xt::range(_hb_range_an.start(), _hb_range_an.end()));
    } else if (_hb_range_ap.valid()) {
      _db = xt::view(gridSpacePtr()->eSizeY(),
                     xt::range(_hb_range_ap.start(), _hb_range_ap.end()));
    }

    if (Axis::directionPositive(_direction)) {
      _positive = 1.0;
    } else {
      _positive = -1.0;
    }
  }
}

void CurrentMonitor::update() {
  auto emf{emfPtr()};
  auto t{calculationParamPtr()->timeParam()->currentTimeStep()};
  double integral_x{0};
  double integral_y{0};
  double integral_z{0};

  if (!valid()) {
    return;
  }

  if (Axis::formDirectionToXYZ(_direction) == Axis::XYZ::X) {
    // for (size_t j{_js}; j <= _je; ++j) {
    //   integral_y +=
    //       (emf->hy()(_ie - 1, j, _ks - 1) - emf->hy()(_ie - 1, j, _ke)) *
    //       _da(j - _js);
    // }
    // for (size_t k{_ks}; k <= _ke; ++k) {
    //   integral_z +=
    //       (emf->hz()(_ie - 1, _je, k) - emf->hz()(_ie - 1, _js - 1, k)) *
    //       _db(k - _ks);
    // }
    // data()(1, t) = _positive * (integral_z + integral_y);
    if (_ha_range_bn.valid()) {
      for (auto j{_ha_range_bn.start()}; j < _ha_range_bn.end(); ++j) {
        integral_y +=
            emf->hy()(_ie - 1, j, _ks - 1) * _da(j - _ha_range_bn.start());
      }
    }
    if (_ha_range_bp.valid()) {
      for (auto j{_ha_range_bp.start()}; j < _ha_range_bp.end(); ++j) {
        integral_y -=
            emf->hy()(_ie - 1, j, _ke) * _da(j - _ha_range_bp.start());
      }
    }
    if (_hb_range_an.valid()) {
      for (auto k{_hb_range_an.start()}; k < _hb_range_an.end(); ++k) {
        integral_z -=
            emf->hz()(_ie - 1, _js - 1, k) * _db(k - _hb_range_an.start());
      }
    }
    if (_hb_range_ap.valid()) {
      for (auto k{_hb_range_ap.start()}; k < _hb_range_ap.end(); ++k) {
        integral_z +=
            emf->hz()(_ie - 1, _je, k) * _db(k - _hb_range_ap.start());
      }
    }
    _node_data(t) = _positive * (integral_z + integral_y);
    return;
  }

  if (Axis::formDirectionToXYZ(_direction) == Axis::XYZ::Y) {
    // for (size_t k{_ks}; k <= _ke; ++k) {
    //   integral_z +=
    //       (emf->hz()(_is - 1, _je - 1, k) - emf->hz()(_ie, _je - 1, k)) *
    //       _da(k - _ks);
    // }
    // for (size_t i{_is}; i <= _ie; ++i) {
    //   integral_x +=
    //       (emf->hx()(i, _je - 1, _ke) - emf->hx()(i, _je - 1, _ks - 1)) *
    //       _db(i - _is);
    // }
    // data()(1, t) = _positive * (integral_x + integral_z);
    if (_ha_range_bn.valid()) {
      for (auto k{_ha_range_bn.start()}; k < _ha_range_bn.end(); ++k) {
        integral_z +=
            emf->hz()(_is - 1, _je - 1, k) * _da(k - _ha_range_bn.start());
      }
    }

    if (_ha_range_bp.valid()) {
      for (auto k{_ha_range_bp.start()}; k < _ha_range_bp.end(); ++k) {
        integral_z -=
            emf->hz()(_ie, _je - 1, k) * _da(k - _ha_range_bp.start());
      }
    }

    if (_hb_range_an.valid()) {
      for (auto i{_hb_range_an.start()}; i < _hb_range_an.end(); ++i) {
        integral_x -=
            emf->hx()(i, _je - 1, _ks - 1) * _db(i - _hb_range_an.start());
      }
    }

    if (_hb_range_ap.valid()) {
      for (auto i{_hb_range_ap.start()}; i < _hb_range_ap.end(); ++i) {
        integral_x +=
            emf->hx()(i, _je - 1, _ke) * _db(i - _hb_range_ap.start());
      }
    }
    _node_data(t) = _positive * (integral_x + integral_z);
    return;
  }

  if (Axis::formDirectionToXYZ(_direction) == Axis::XYZ::Z) {
    if (_ha_range_bn.valid()) {
      for (auto i{_ha_range_bn.start()}; i < _ha_range_bn.end(); ++i) {
        integral_x +=
            (emf->hx().at(i, _js - 1, _ke - 1)) * _da(i - _ha_range_bn.start());
      }
    }

    if (_ha_range_bp.valid()) {
      for (auto i{_ha_range_bp.start()}; i < _ha_range_bp.end(); ++i) {
        integral_x -=
            (emf->hx()(i, _je, _ke - 1)) * _da(i - _ha_range_bp.start());
      }
    }

    if (_hb_range_an.valid()) {
      for (auto j{_hb_range_an.start()}; j < _hb_range_an.end(); ++j) {
        integral_y -=
            (emf->hy()(_is - 1, j, _ke - 1)) * _db(j - _hb_range_an.start());
      }
    }

    if (_hb_range_ap.valid()) {
      for (auto j{_hb_range_ap.start()}; j < _hb_range_ap.end(); ++j) {
        integral_y +=
            (emf->hy()(_ie, j, _ke - 1)) * _db(j - _hb_range_ap.start());
      }
    }

    _node_data(t) = _positive * (integral_x + integral_y);
    return;
  }
}

void CurrentMonitor::initTimeDependentVariable() {
  setTime(calculationParamPtr()->timeParam()->hTime());
  _node_data = xt::empty_like(time());
  data() = xt::empty_like(time());
}

void CurrentMonitor::initParallelizedConfig() { makeMpiSubComm(); }

auto CurrentMonitor::gatherData() -> void {
  if (!valid()) {
    return;
  }

  if (monitorMpiConfig().size() <= 1) {
    data() = _node_data;
    return;
  }

  data() = xt::empty_like(time());

  auto& mpi_support = MpiSupport::instance();

  mpi_support.reduceSum(monitorMpiConfig(), _node_data.data(), data().data(),
                        data().size());
}

auto CurrentMonitor::valid() const -> bool {
  return _ha_range_bn.valid() || _ha_range_bp.valid() || _hb_range_an.valid() ||
         _hb_range_ap.valid();
}

auto CurrentMonitor::toString() const -> std::string {
  std::stringstream ss;
  ss << "Current Monitor:\n";
  ss << " Global " << globalGridBox().toString() << "\n";
  ss << " Node " << nodeGridBox().toString() << "\n";
  // ss << " Direction " << Axis::directionToString(_direction) << "\n";
  ss << " Ha Range Bn " << _ha_range_bn.toString() << "\n";
  ss << " Ha Range Bp " << _ha_range_bp.toString() << "\n";
  ss << " Hb Range An " << _hb_range_an.toString() << "\n";
  ss << " Hb Range Ap " << _hb_range_ap.toString() << "\n";
  return ss.str();
}

}  // namespace xfdtd
