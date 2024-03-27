#include <xfdtd/parallel/mpi_support.h>

#include <sstream>
#include <string>
#include <vector>

#if defined(XFDTD_CORE_WITH_MPI)
#include <mpi.h>
#endif

namespace xfdtd {

MpiSupport::MpiSupport(int argc, char** argv) { mpiInit(argc, argv); }

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

#if defined(XFDTD_CORE_WITH_MPI)

// auto MpiSupport::barrier(MPI_Comm comm) const -> void { MPI_Barrier(comm); }
#endif

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

auto MpiSupport::sendRecvHyXHead(xt::xarray<double>& hy) -> void {
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

auto MpiSupport::recvSendHyXTail(xt::xarray<double>& hy) -> void {
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

auto MpiSupport::sendRecvHzXHead(xt::xarray<double>& hz) -> void {
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

auto MpiSupport::recvSendHzXTail(xt::xarray<double>& hz) -> void {
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

auto MpiSupport::sendRecvHzYHead(xt::xarray<double>& hz) -> void {
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

auto MpiSupport::recvSendHzYTail(xt::xarray<double>& hz) -> void {
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

auto MpiSupport::sendRecvHxYHead(xt::xarray<double>& hx) -> void {
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

auto MpiSupport::recvSendHxYTail(xt::xarray<double>& hx) -> void {
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

auto MpiSupport::sendRecvHxZHead(xt::xarray<double>& hx) -> void {
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

auto MpiSupport::recvSendHxZTail(xt::xarray<double>& hx) -> void {
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

auto MpiSupport::sendRecvHyZHead(xt::xarray<double>& hy) -> void {
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

auto MpiSupport::recvSendHyZTail(xt::xarray<double>& hy) -> void {
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

  int num_x = 1;
  int num_y = 1;
  int num_z = 1;

  num_x = _global_size;
  if (4 <= argc) {
    if (_global_rank == 0) {
      std::stringstream ss;
      ss << "Warning: Mpi Support will read command line arguments. Usage: "
            " mpiexec -n <num_procs> ./<executable> <num_x> <num_y> <num_z>\n";
      ss << "Warning: Here is temporary way to set layout of MPI ranks.\n";
      std::cout << ss.str();
    }

    num_x = std::stoi(argv[1]);
    num_y = std::stoi(argv[2]);
    num_z = std::stoi(argv[3]);
  }

  std::vector<int> ranks(_global_size, 0);
  std::iota(ranks.begin(), ranks.end(), 0);
  _config = MpiConfig::makeXFDTDComm(_global_rank, 0, num_x, num_y, num_z);

#endif
}

MpiSupport::MpiGuard::~MpiGuard() {
#if defined(XFDTD_CORE_WITH_MPI)
  MPI_Finalize();
#endif
}

auto MpiSupport::toString() const -> std::string { return _config.toString(); }

}  // namespace xfdtd
