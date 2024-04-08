#include <xfdtd/common/type_define.h>
#include <xfdtd/network/port.h>
#include <xfdtd/parallel/mpi_support.h>
#include <xfdtd/util/dft.h>

#include <complex>

namespace xfdtd {

Port::Port(std::size_t index, bool is_source, std::complex<Real> impedance,
           std::shared_ptr<CurrentMonitor> current_monitor,
           std::shared_ptr<VoltageMonitor> voltage_monitor)
    : _index{index},
      _is_source{is_source},
      _impedance{impedance},
      _current_monitor{std::move(current_monitor)},
      _voltage_monitor{std::move(voltage_monitor)} {}

std::size_t Port::index() const { return _index; }

bool Port::isSource() const { return _is_source; }

std::complex<Real> Port::impedance() const { return _impedance; }

const std::shared_ptr<CurrentMonitor>& Port::currentMonitor() const {
  return _current_monitor;
}

const std::shared_ptr<VoltageMonitor>& Port::voltageMonitor() const {
  return _voltage_monitor;
}

const Array1D<std::complex<Real>>& Port::a() const { return _a; }

const Array1D<std::complex<Real>>& Port::b() const { return _b; }

void Port::init(
    const std::shared_ptr<const GridSpace>& grid_space,
    const std::shared_ptr<const CalculationParam>& calculation_param,
    const std::shared_ptr<const EMF>& emf) {
  _dt = calculation_param->timeParam()->dt();
}

void Port::calculateSParameters(const Array1D<Real>& frequencies) {
  const auto& c_const = *_current_monitor;
  const auto& v_const = *_voltage_monitor;

  Array1D<Real> current = xt::zeros_like(c_const.data());
  Array1D<Real> voltage = xt::zeros_like(v_const.data());

  /**
   * @brief Only root in the xfdtd comm can write data to file. First, let the
   * root in the monitor comm gather data, then send it to the root in the xfdtd
   * comm. Communicator is xfdtd and using tag to distinguish current and
   * voltage. Careful: xfdtd root don't know who is root in the monitor comm, so
   * we use MPI_ANY_SOURCE to recv mes.
   *
   */
  auto collect_func = [](const auto& const_monitor, auto& monitor, auto& data,
                         const auto& tag) {
    auto& mpi_support = MpiSupport::instance();

    if (!const_monitor.valid()) {
      if (!mpi_support.isRoot()) {
        return;
      }

      if (mpi_support.size() <= 1) {
        throw XFDTDException(
            "Port: Monitor is not valid, but MPI size is not 1. Please make "
            "sure you had add monitor to the simulation.");
      }

      // xfdtd root is not in the monitor comm.
      mpi_support.recv(mpi_support.config(), data.data(),
                       sizeof(Real) * data.size(), MpiSupport::ANY_SOURCE, tag);
      return;
    }

    monitor.gatherData();
    if (const_monitor.monitorMpiConfig().isRoot()) {
      data = monitor.data();
    }

    if (mpi_support.size() <= 1) {
      return;
    }

    if (!mpi_support.isRoot() && const_monitor.monitorMpiConfig().isRoot()) {
      // monitor root send data to xfdtd root.
      mpi_support.send(mpi_support.config(), data.data(),
                       sizeof(Real) * data.size(), mpi_support.config().root(),
                       tag);
      return;
    }

    if (mpi_support.isRoot() && !const_monitor.monitorMpiConfig().isRoot()) {
      // xfdtd root is in the monitor comm, but not the root of the monitor
      // comm.
      mpi_support.recv(mpi_support.config(), data.data(),
                       sizeof(Real) * data.size(), MpiSupport::ANY_SOURCE, tag);
    }

    // 1. not root in the xfdtd comm and not root in the monitor comm.
    // 2. root in the xfdtd comm and root in the monitor comm.
    // These two cases do not need to do anything.
  };

  // multi-threading or single-threading
  // std::vector<std::thread> collector;

  collect_func(c_const, *_current_monitor, current, CURRENT_TAG);
  collect_func(v_const, *_voltage_monitor, voltage, VOLTAGE_TAG);

  MpiSupport::instance().barrier();

  if (!MpiSupport::instance().isRoot()) {
    return;
  }

  auto dt{_dt};
  auto z{_impedance};
  auto k{std::sqrt(std::real(z))};
  auto i{dft(current, dt, frequencies, -0.5 * dt)};
  auto v{dft(voltage, dt, frequencies)};

  _a = 0.5 * (v + z * i) / k;
  _b = 0.5 * (v - std::conj(z) * i) / k;
}

}  // namespace xfdtd
