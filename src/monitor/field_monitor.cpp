#include <xfdtd/grid_space/grid_space.h>
#include <xfdtd/monitor/field_monitor.h>
#include <xfdtd/monitor/monitor.h>

#include <utility>
#include <xtensor.hpp>

#include "xfdtd/parallel/mpi_support.h"

namespace xfdtd {

FieldMonitor::FieldMonitor(std::unique_ptr<Shape> shape, EMF::Field field,
                           std::string name, std::string output_dir_path)
    : Monitor{std::move(shape), std::move(name), std::move(output_dir_path)},
      _field{field} {}

void FieldMonitor::init(
    std::shared_ptr<const GridSpace> grid_space,
    std::shared_ptr<const CalculationParam> calculation_param,
    std::shared_ptr<const EMF> emf) {
  defaultInit(grid_space, calculation_param, emf);

  if (gridSpacePtr()->dimension() == GridSpace::Dimension::ONE) {
    throw XFDTDMonitorException(
        "FieldMonitor cannot be used in 1D simulation(not implemented)");
  }

  // TODO(franzero): temporary way. need to be refactored
  // Example: The box for Hx is. The box for Ex is
  auto offset_i = 0;
  auto offset_j = 0;
  auto offset_k = 0;

  if (globalGridBox().origin().i() == 0) {
    offset_i = 1;
  }
  if (globalGridBox().origin().j() == 0) {
    offset_j = 1;
  }
  if (globalGridBox().origin().k() == 0) {
    offset_k =
        gridSpacePtr()->dimension() == GridSpace::Dimension::THREE ? 1 : 0;
  }
  setGlobalGridBox(
      GridBox{globalGridBox().origin() + Grid{static_cast<size_t>(offset_i),
                                              static_cast<size_t>(offset_j),
                                              static_cast<size_t>(offset_k)},
              globalGridBox().size() - Grid{static_cast<size_t>(offset_i),
                                            static_cast<size_t>(offset_j),
                                            static_cast<size_t>(offset_k)}});
  offset_i = 0;
  offset_j = 0;
  offset_k = 0;

  if (nodeGridBox().origin().i() == 0) {
    offset_i = 1;
  }
  if (nodeGridBox().origin().j() == 0) {
    offset_j = 1;
  }
  if (nodeGridBox().origin().k() == 0) {
    offset_k =
        gridSpacePtr()->dimension() == GridSpace::Dimension::THREE ? 1 : 0;
  }

  setNodeGridBox(
      GridBox{nodeGridBox().origin() + Grid{static_cast<size_t>(offset_i),
                                            static_cast<size_t>(offset_j),
                                            static_cast<size_t>(offset_k)},
              nodeGridBox().size() - Grid{static_cast<size_t>(offset_i),
                                          static_cast<size_t>(offset_j),
                                          static_cast<size_t>(offset_k)}});

  setGlobalTask(makeIndexTask(
      makeIndexRange(globalGridBox().origin().i(), globalGridBox().end().i()),
      makeIndexRange(globalGridBox().origin().j(), globalGridBox().end().j()),
      makeIndexRange(globalGridBox().origin().k(), globalGridBox().end().k())));

  setNodeTask(makeIndexTask(
      makeIndexRange(nodeGridBox().origin().i(), nodeGridBox().end().i()),
      makeIndexRange(nodeGridBox().origin().j(), nodeGridBox().end().j()),
      makeIndexRange(nodeGridBox().origin().k(), nodeGridBox().end().k())));
}

void FieldMonitor::update() {}

void FieldMonitor::output() {
  if (!valid()) {
    return;
  }
  auto em_field{emfPtr()};

  auto grid_box{nodeGridBox()};
  auto grid_box_size{grid_box.size()};
  auto x_range{xt::range(grid_box.origin().i(), grid_box.end().i())};
  auto y_range{xt::range(grid_box.origin().j(), grid_box.end().j())};
  auto z_range{xt::range(grid_box.origin().k(), grid_box.end().k())};
  data() = xt::view(em_field->field(field()), x_range, y_range, z_range);
  gatherData();
  Monitor::output();
}

EMF::Field FieldMonitor::field() const { return _field; }

auto FieldMonitor::initParallelizedConfig() -> void {
  makeMpiSubComm();
  if (!valid()) {
    return;
  }

  auto& mpi_support = MpiSupport::instance();

#if defined(XFDTD_CORE_WITH_MPI)
  const auto node_box_origin_in_global =
      gridSpacePtr()->globalBox().origin() + nodeGridBox().origin();
  const auto g_origin = globalGridBox().origin();

  auto nx = nodeTask().xRange().size();
  auto ny = nodeTask().yRange().size();
  auto nz = nodeTask().zRange().size();
  auto stride_elem = globalTask().zRange().size();
  auto stride_vec = globalTask().zRange().size() * globalTask().yRange().size();
  auto disp = (node_box_origin_in_global.i() - g_origin.i()) * stride_vec +
              (node_box_origin_in_global.j() - g_origin.j()) * stride_elem +
              (node_box_origin_in_global.k() - g_origin.k());

  auto p = MpiSupport::Block::Profile{
      static_cast<int>(nx),          static_cast<int>(ny),
      static_cast<int>(nz),          static_cast<int>(stride_vec),
      static_cast<int>(stride_elem), static_cast<int>(disp)};

  if (monitorMpiConfig().isRoot()) {
    _profiles.resize(monitorMpiConfig().size());
  }

  auto profile_type = MpiSupport::TypeGuard{};
  MPI_Type_contiguous(sizeof(MpiSupport::Block::Profile), MPI_CHAR,
                      &profile_type._type);
  MPI_Type_commit(&profile_type._type);

  mpi_support.gather(monitorMpiConfig(), &p, 1, profile_type, _profiles.data(),
                     1, profile_type, monitorMpiConfig().root());

  if (monitorMpiConfig().isRoot()) {
    _blocks_mpi.reserve(_profiles.size());
    for (const auto& profile : _profiles) {
      _blocks_mpi.emplace_back(MpiSupport::Block::make(profile));
    }
  }
#endif
}

auto FieldMonitor::gatherData() -> void {
  if (!valid() || monitorMpiConfig().size() == 1) {
    return;
  }

  auto& mpi_support = MpiSupport::instance();
  if (monitorMpiConfig().isRoot()) {
    Array3D<Real> recv_buffer = xt::zeros<Real>({globalTask().xRange().size(),
                                                 globalTask().yRange().size(),
                                                 globalTask().zRange().size()});
    for (int i = 1; i < monitorMpiConfig().size(); ++i) {
      mpi_support.iRecv(monitorMpiConfig(), recv_buffer.data(), 1,
                        _blocks_mpi[i], i, 0);
    }

    mpi_support.waitAll();

    mpi_support.sendRecv(monitorMpiConfig(), data().data(), data().size(),
                         monitorMpiConfig().rank(), 0, recv_buffer.data(), 1,
                         _blocks_mpi[0], monitorMpiConfig().rank(), 0);

    data() = std::move(recv_buffer);
  } else {
    mpi_support.send(monitorMpiConfig(), data().data(),
                     sizeof(Real) * data().size(), monitorMpiConfig().root(),
                     0);
  }
}

}  // namespace xfdtd
