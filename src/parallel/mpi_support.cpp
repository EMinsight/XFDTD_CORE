#include <xfdtd/parallel/mpi_support.h>

#include <sstream>
#include <string>
#include <vector>

#if defined(XFDTD_CORE_WITH_MPI)
#include <mpi.h>
#endif

namespace xfdtd {

auto MpiSupport::setMpiParallelDim(int nx, int ny, int nz) -> void {
  static bool is_set = false;
  if (is_set) {
    return;
  }
  config_nx = nx;
  config_ny = ny;
  config_nz = nz;
}

MpiSupport::MpiSupport(int argc, char** argv) {
#if defined(XFDTD_CORE_SINGLE_PRECISION)
  throw XFDTDException("MpiSupport is not supported with single precision.");
#endif
  mpiInit(argc, argv);
}

auto MpiSupport::abort(int error_code) const -> void {
#if defined(XFDTD_CORE_WITH_MPI)
  MPI_Abort(_config.comm(), error_code);
#else
  std::exit(error_code);
#endif
}

auto MpiSupport::barrier() const -> void {
#if defined(XFDTD_CORE_WITH_MPI)
  MPI_Barrier(_config.comm());
#endif
}

auto MpiSupport::generateSlice(std::size_t nx, std::size_t ny, std::size_t nz)
    -> void {
  _hx_z_slice = Slice::makeHxZSlice(nx, ny, nz);
  _hy_z_slice = Slice::makeHyZSlice(nx, ny, nz);

  _hx_y_slice = Slice::makeHxYSlice(nx, ny, nz);
  _hz_y_slice = Slice::makeHzYSlice(nx, ny, nz);

  _hy_x_slice = Slice::makeHyXSlice(nx, ny, nz);
  _hz_x_slice = Slice::makeHzXSlice(nx, ny, nz);
}

auto MpiSupport::xNext() const -> int { return _config.xNext(); }

auto MpiSupport::xPrev() const -> int { return _config.xPrev(); }

auto MpiSupport::yNext() const -> int { return _config.yNext(); }

auto MpiSupport::yPrev() const -> int { return _config.yPrev(); }

auto MpiSupport::zNext() const -> int { return _config.zNext(); }

auto MpiSupport::zPrev() const -> int { return _config.zPrev(); }

auto MpiSupport::sendRecvHyXHead(Array3D<Real>& hy) -> void {
#if defined(XFDTD_CORE_WITH_MPI)
  MPI_Request s_request;
  MPI_Request r_request;

  MPI_Isend(&hy(1, 0, 0), 1, _hy_x_slice.slice(), xPrev(), exchange_hy_x_sr_tag,
            _config.comm(), &s_request);
  MPI_Irecv(hy.data(), 1, _hy_x_slice.slice(), xPrev(), exchange_hy_x_rs_tag,
            _config.comm(), &r_request);

  _requests.emplace_back(s_request);
  _requests.emplace_back(r_request);
#endif
}

auto MpiSupport::recvSendHyXTail(Array3D<Real>& hy) -> void {
#if defined(XFDTD_CORE_WITH_MPI)
  const auto nx = hy.shape()[0];

  MPI_Request r_request;
  MPI_Request s_request;

  MPI_Irecv(&hy(nx - 1, 0, 0), 1, _hy_x_slice.slice(), xNext(),
            exchange_hy_x_sr_tag, _config.comm(), &r_request);
  MPI_Isend(&hy(nx - 2, 0, 0), 1, _hy_x_slice.slice(), xNext(),
            exchange_hy_x_rs_tag, _config.comm(), &s_request);

  _requests.emplace_back(r_request);
  _requests.emplace_back(s_request);
#endif
}

auto MpiSupport::sendRecvHzXHead(Array3D<Real>& hz) -> void {
#if defined(XFDTD_CORE_WITH_MPI)
  MPI_Request s_request;
  MPI_Request r_request;

  MPI_Isend(&hz(1, 0, 0), 1, _hz_x_slice.slice(), xPrev(), exchange_hz_x_sr_tag,
            _config.comm(), &s_request);
  MPI_Irecv(hz.data(), 1, _hz_x_slice.slice(), xPrev(), exchange_hz_x_rs_tag,
            _config.comm(), &r_request);

  _requests.emplace_back(s_request);
  _requests.emplace_back(r_request);
#endif
}

auto MpiSupport::recvSendHzXTail(Array3D<Real>& hz) -> void {
#if defined(XFDTD_CORE_WITH_MPI)
  const auto nx = hz.shape()[0];

  MPI_Request r_request;
  MPI_Request s_request;

  MPI_Irecv(&hz(nx - 1, 0, 0), 1, _hz_x_slice.slice(), xNext(),
            exchange_hz_x_sr_tag, _config.comm(), &r_request);
  MPI_Isend(&hz(nx - 2, 0, 0), 1, _hz_x_slice.slice(), xNext(),
            exchange_hz_x_rs_tag, _config.comm(), &s_request);

  _requests.emplace_back(r_request);
  _requests.emplace_back(s_request);
#endif
}

auto MpiSupport::sendRecvHzYHead(Array3D<Real>& hz) -> void {
#if defined(XFDTD_CORE_WITH_MPI)
  MPI_Request s_request;
  MPI_Request r_request;

  MPI_Isend(&hz(0, 1, 0), 1, _hz_y_slice.slice(), yPrev(), exchange_hz_y_sr_tag,
            _config.comm(), &s_request);
  MPI_Irecv(hz.data(), 1, _hz_y_slice.slice(), yPrev(), exchange_hz_y_rs_tag,
            _config.comm(), &r_request);

  _requests.emplace_back(s_request);
  _requests.emplace_back(r_request);
#endif
}

auto MpiSupport::recvSendHzYTail(Array3D<Real>& hz) -> void {
#if defined(XFDTD_CORE_WITH_MPI)
  const auto ny = hz.shape()[1];

  MPI_Request r_request;
  MPI_Request s_request;

  MPI_Irecv(&hz(0, ny - 1, 0), 1, _hz_y_slice.slice(), yNext(),
            exchange_hz_y_sr_tag, _config.comm(), &r_request);
  MPI_Isend(&hz(0, ny - 2, 0), 1, _hz_y_slice.slice(), yNext(),
            exchange_hz_y_rs_tag, _config.comm(), &s_request);

  _requests.emplace_back(r_request);
  _requests.emplace_back(s_request);
#endif
}

auto MpiSupport::sendRecvHxYHead(Array3D<Real>& hx) -> void {
#if defined(XFDTD_CORE_WITH_MPI)
  MPI_Request s_request;
  MPI_Request r_request;

  MPI_Isend(&hx(0, 1, 0), 1, _hx_y_slice.slice(), yPrev(), exchange_hx_y_sr_tag,
            _config.comm(), &s_request);
  MPI_Irecv(hx.data(), 1, _hx_y_slice.slice(), yPrev(), exchange_hx_y_rs_tag,
            _config.comm(), &r_request);

  _requests.emplace_back(s_request);
  _requests.emplace_back(r_request);
#endif
}

auto MpiSupport::recvSendHxYTail(Array3D<Real>& hx) -> void {
#if defined(XFDTD_CORE_WITH_MPI)
  const auto ny = hx.shape()[1];

  MPI_Request r_request;
  MPI_Request s_request;

  MPI_Irecv(&hx(0, ny - 1, 0), 1, _hx_y_slice.slice(), yNext(),
            exchange_hx_y_sr_tag, _config.comm(), &r_request);
  MPI_Isend(&hx(0, ny - 2, 0), 1, _hx_y_slice.slice(), yNext(),
            exchange_hx_y_rs_tag, _config.comm(), &s_request);

  _requests.emplace_back(r_request);
  _requests.emplace_back(s_request);
#endif
}

auto MpiSupport::sendRecvHxZHead(Array3D<Real>& hx) -> void {
#if defined(XFDTD_CORE_WITH_MPI)
  MPI_Request s_request;
  MPI_Request r_request;

  MPI_Isend(&hx(0, 0, 1), 1, _hx_z_slice.slice(), zPrev(), exchange_hx_z_sr_tag,
            _config.comm(), &s_request);
  MPI_Irecv(hx.data(), 1, _hx_z_slice.slice(), zPrev(), exchange_hx_z_rs_tag,
            _config.comm(), &r_request);

  _requests.emplace_back(s_request);
  _requests.emplace_back(r_request);
#endif
}

auto MpiSupport::recvSendHxZTail(Array3D<Real>& hx) -> void {
#if defined(XFDTD_CORE_WITH_MPI)
  const auto nz = hx.shape()[2];

  MPI_Request r_request;
  MPI_Request s_request;

  MPI_Irecv(&hx(0, 0, nz - 1), 1, _hx_z_slice.slice(), zNext(),
            exchange_hx_z_sr_tag, _config.comm(), &r_request);
  MPI_Isend(&hx(0, 0, nz - 2), 1, _hx_z_slice.slice(), zNext(),
            exchange_hx_z_rs_tag, _config.comm(), &s_request);

  _requests.emplace_back(r_request);
  _requests.emplace_back(s_request);
#endif
}

auto MpiSupport::sendRecvHyZHead(Array3D<Real>& hy) -> void {
#if defined(XFDTD_CORE_WITH_MPI)
  MPI_Request s_request;
  MPI_Request r_request;

  MPI_Isend(&hy(0, 0, 1), 1, _hy_z_slice.slice(), zPrev(), exchange_hy_z_sr_tag,
            _config.comm(), &s_request);
  MPI_Irecv(hy.data(), 1, _hy_z_slice.slice(), zPrev(), exchange_hy_z_rs_tag,
            _config.comm(), &r_request);

  _requests.emplace_back(s_request);
  _requests.emplace_back(r_request);
#endif
}

auto MpiSupport::recvSendHyZTail(Array3D<Real>& hy) -> void {
#if defined(XFDTD_CORE_WITH_MPI)
  const auto nz = hy.shape()[2];

  MPI_Request r_request;
  MPI_Request s_request;

  MPI_Irecv(&hy(0, 0, nz - 1), 1, _hy_z_slice.slice(), zNext(),
            exchange_hy_z_sr_tag, _config.comm(), &r_request);
  MPI_Isend(&hy(0, 0, nz - 2), 1, _hy_z_slice.slice(), zNext(),
            exchange_hy_z_rs_tag, _config.comm(), &s_request);

  _requests.emplace_back(r_request);
  _requests.emplace_back(s_request);
#endif
}

auto MpiSupport::waitAll() -> void {
#if defined(XFDTD_CORE_WITH_MPI)
  std::vector<MPI_Status> statuses(_requests.size());
  MPI_Waitall(_requests.size(), _requests.data(), statuses.data());
  _requests.clear();
#endif
}

auto MpiSupport::mpiInit(int argc, char** argv) -> void {
#if defined(XFDTD_CORE_WITH_MPI)
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &_global_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &_global_size);

  if (_global_size != config_nx * config_ny * config_nz) {
    if (_global_rank == 0) {
      std::stringstream ss;
      ss << "Info: MPI configuration "
         << "(" << config_nx << " x " << config_nx << " x " << config_nx << ") "
         << "is not matched with the number of "
            "processes. "
            "Set the layout of MPI configuration to "
         << _global_size << " x " << 1 << " x " << 1 << ".\n";

      std::cout << ss.str();
    }
    config_nx = _global_size;
    config_ny = 1;
    config_nz = 1;
  }

  _config = MpiConfig::makeXFDTDComm(_global_rank, 0, config_nx, config_ny,
                                     config_nz);

  if (isRoot()) {
    std::stringstream ss;
    ss << "Info: MPI configuration "
       << "(" << config_nx << " x " << config_nx << " x " << config_nx << ") "
       << "is set.\n";
    std::cout << ss.str();
  }

#endif
}

MpiSupport::MpiGuard::~MpiGuard() {
#if defined(XFDTD_CORE_WITH_MPI)
  MPI_Finalize();
#endif
}

auto MpiSupport::toString() const -> std::string { return _config.toString(); }

}  // namespace xfdtd
