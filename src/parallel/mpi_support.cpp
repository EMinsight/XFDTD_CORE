#include <xfdtd/electromagnetic_field/electromagnetic_field.h>
#include <xfdtd/parallel/mpi_support.h>

#include <cstddef>
#include <cstdlib>
#include <sstream>
#include <string>
#include <vector>

#include "xfdtd/common/type_define.h"

#if defined(XFDTD_CORE_WITH_MPI)
#include <mpi.h>
#endif

namespace xfdtd {

MpiSupport::~MpiSupport() { waitAll(); }

auto MpiSupport::setMpiParallelDim(int nx, int ny, int nz) -> void {
  static bool is_set = false;
  if (is_set) {
    return;
  }
  config_nx = nx;
  config_ny = ny;
  config_nz = nz;
}

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

inline static auto getAddress(const Real* ptr) -> ptrdiff_t {
  return reinterpret_cast<ptrdiff_t>(ptr);
}

inline static auto getRealAddressOffset(const Real* ptr_a, const Real* ptr_b)
    -> ptrdiff_t {
  return getAddress(ptr_a) - getAddress(ptr_b);
}

auto MpiSupport::generateSlice(std::size_t nx, std::size_t ny, std::size_t nz)
    -> void {
  _hy_xn_s_block = Block::makeRowMajorXSlice(1, ny + 1, nz, (ny + 1) * nz);
  _hy_xn_r_block = Block::makeRowMajorXSlice(1, ny + 1, nz, 0);

  _hy_xp_r_block =
      Block::makeRowMajorXSlice(1, ny + 1, nz, ((ny + 1) * nz) * (nx - 1));
  _hy_xp_s_block =
      Block::makeRowMajorXSlice(1, ny + 1, nz, ((ny + 1) * nz) * (nx - 2));

  _hz_xn_s_block = Block::makeRowMajorXSlice(1, ny, nz + 1, ny * (nz + 1));
  _hz_xn_r_block = Block::makeRowMajorXSlice(1, ny, nz + 1, 0);

  _hz_xp_r_block =
      Block::makeRowMajorXSlice(1, ny, nz + 1, (ny * (nz + 1)) * (nx - 1));
  _hz_xp_s_block =
      Block::makeRowMajorXSlice(1, ny, nz + 1, (ny * (nz + 1)) * (nx - 2));

  _hz_yn_s_block =
      Block::makeRowMajorYSlice(1, nx, ny, nz + 1, (nz + 1));  // send hz(0,1,0)
  _hz_yn_r_block =
      Block::makeRowMajorYSlice(1, nx, ny, nz + 1, 0);  // recv hz(0,0,0)

  _hz_yp_r_block = Block::makeRowMajorYSlice(
      1, nx, ny, nz + 1, (nz + 1) * (ny - 1));  // recv hz(0,ny-1,0)
  _hz_yp_s_block = Block::makeRowMajorYSlice(
      1, nx, ny, nz + 1, (nz + 1) * (ny - 2));  // send hz(0,ny-2,0)

  _hx_yn_s_block =
      Block::makeRowMajorYSlice(1, nx + 1, ny, nz, nz);  // send hx(0,1,0)
  _hx_yn_r_block =
      Block::makeRowMajorYSlice(1, nx + 1, ny, nz, 0);  // recv hx(0,0,0)

  _hx_yp_r_block =
      Block::makeRowMajorYSlice(1, nx + 1, ny, nz,
                                (ny - 1) * nz);  // recv hx(0,ny-1,0)
  _hx_yp_s_block =
      Block::makeRowMajorYSlice(1, nx + 1, ny, nz,
                                (ny - 2) * nz);  // send hx(0,ny-2,0)

  _hx_zn_s_block =
      Block::makeRowMajorZSlice(1, nx + 1, ny, nz, 1);  // send hx(0,0,1)
  _hx_zn_r_block =
      Block::makeRowMajorZSlice(1, nx + 1, ny, nz, 0);  // recv hx(0,0,0)

  _hx_zp_r_block = Block::makeRowMajorZSlice(1, nx + 1, ny, nz,
                                             (nz - 1));  // recv hx(0,0,nz-1)
  _hx_zp_s_block = Block::makeRowMajorZSlice(1, nx + 1, ny, nz,
                                             (nz - 2));  // send hx(0,0,nz-2)

  _hy_zn_s_block =
      Block::makeRowMajorZSlice(1, nx, ny + 1, nz, 1);  // send hy(0,0,1)
  _hy_zn_r_block =
      Block::makeRowMajorZSlice(1, nx, ny + 1, nz, 0);  // recv hy(0,0,0)

  _hy_zp_r_block = Block::makeRowMajorZSlice(1, nx, ny + 1, nz,
                                             (nz - 1));  // recv hy(0,0,nz-1)
  _hy_zp_s_block = Block::makeRowMajorZSlice(1, nx, ny + 1, nz,
                                             (nz - 2));  // send hy(0,0,nz-2)
}

auto MpiSupport::xNext() const -> int { return _config.xNext(); }

auto MpiSupport::xPrev() const -> int { return _config.xPrev(); }

auto MpiSupport::yNext() const -> int { return _config.yNext(); }

auto MpiSupport::yPrev() const -> int { return _config.yPrev(); }

auto MpiSupport::zNext() const -> int { return _config.zNext(); }

auto MpiSupport::zPrev() const -> int { return _config.zPrev(); }

auto MpiSupport::sendRecvHyXHead(Array3D<Real>& hy) -> void {
  iSend(config(), hy.data(), 1, _hy_xn_s_block, xPrev(), EXCHANGE_HY_X_SR_TAG);
  iRecv(config(), hy.data(), 1, _hy_xn_r_block, xPrev(), EXCHANGE_HY_X_RS_TAG);
}

auto MpiSupport::recvSendHyXTail(Array3D<Real>& hy) -> void {
  iRecv(config(), hy.data(), 1, _hy_xp_r_block, xNext(), EXCHANGE_HY_X_SR_TAG);
  iSend(config(), hy.data(), 1, _hy_xp_s_block, xNext(), EXCHANGE_HY_X_RS_TAG);
}

auto MpiSupport::sendRecvHzXHead(Array3D<Real>& hz) -> void {
  iSend(config(), hz.data(), 1, _hz_xn_s_block, xPrev(), EXCHANGE_HZ_X_SR_TAG);
  iRecv(config(), hz.data(), 1, _hz_xn_r_block, xPrev(), EXCHANGE_HZ_X_RS_TAG);
}

auto MpiSupport::recvSendHzXTail(Array3D<Real>& hz) -> void {
  iRecv(config(), hz.data(), 1, _hz_xp_r_block, xNext(), EXCHANGE_HZ_X_SR_TAG);
  iSend(config(), hz.data(), 1, _hz_xp_s_block, xNext(), EXCHANGE_HZ_X_RS_TAG);
}

auto MpiSupport::sendRecvHzYHead(Array3D<Real>& hz) -> void {
  iSend(config(), hz.data(), 1, _hz_yn_s_block, yPrev(), EXCHANGE_HZ_Y_SR_TAG);
  iRecv(config(), hz.data(), 1, _hz_yn_r_block, yPrev(), EXCHANGE_HZ_Y_RS_TAG);
}

auto MpiSupport::recvSendHzYTail(Array3D<Real>& hz) -> void {
  iRecv(config(), hz.data(), 1, _hz_yp_r_block, yNext(), EXCHANGE_HZ_Y_SR_TAG);
  iSend(config(), hz.data(), 1, _hz_yp_s_block, yNext(), EXCHANGE_HZ_Y_RS_TAG);
}

auto MpiSupport::sendRecvHxYHead(Array3D<Real>& hx) -> void {
  iSend(config(), hx.data(), 1, _hx_yn_s_block, yPrev(), EXCHANGE_HX_Y_SR_TAG);
  iRecv(config(), hx.data(), 1, _hx_yn_r_block, yPrev(), EXCHANGE_HX_Y_RS_TAG);
}

auto MpiSupport::recvSendHxYTail(Array3D<Real>& hx) -> void {
  iRecv(config(), hx.data(), 1, _hx_yp_r_block, yNext(), EXCHANGE_HX_Y_SR_TAG);
  iSend(config(), hx.data(), 1, _hx_yp_s_block, yNext(), EXCHANGE_HX_Y_RS_TAG);
}

auto MpiSupport::sendRecvHxZHead(Array3D<Real>& hx) -> void {
  iSend(config(), hx.data(), 1, _hx_zn_s_block, zPrev(), EXCHANGE_HX_Z_SR_TAG);
  iRecv(config(), hx.data(), 1, _hx_zn_r_block, zPrev(), EXCHANGE_HX_Z_RS_TAG);
}

auto MpiSupport::recvSendHxZTail(Array3D<Real>& hx) -> void {
  iRecv(config(), hx.data(), 1, _hx_zp_r_block, zNext(), EXCHANGE_HX_Z_SR_TAG);
  iSend(config(), hx.data(), 1, _hx_zp_s_block, zNext(), EXCHANGE_HX_Z_RS_TAG);
}

auto MpiSupport::sendRecvHyZHead(Array3D<Real>& hy) -> void {
  iSend(config(), hy.data(), 1, _hy_zn_s_block, zPrev(), EXCHANGE_HY_Z_SR_TAG);
  iRecv(config(), hy.data(), 1, _hy_zn_r_block, zPrev(), EXCHANGE_HY_Z_RS_TAG);
}

auto MpiSupport::recvSendHyZTail(Array3D<Real>& hy) -> void {
  iRecv(config(), hy.data(), 1, _hy_zp_r_block, zNext(), EXCHANGE_HY_Z_SR_TAG);
  iSend(config(), hy.data(), 1, _hy_zp_s_block, zNext(), EXCHANGE_HY_Z_RS_TAG);
}

auto MpiSupport::waitAll() -> void {
#if defined(XFDTD_CORE_WITH_MPI)
  std::vector<MPI_Status> statuses(_requests.size());
  MPI_Waitall(_requests.size(), _requests.data(), statuses.data());
  for (const auto& s : statuses) {
    if (s.MPI_ERROR != MPI_SUCCESS) {
      std::cerr << "MPI_Waitall failed\n";
      abort(s.MPI_ERROR);
    }
  }
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
         << "(" << config_nx << " x " << config_ny << " x " << config_nz << ") "
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
       << "(" << config_nx << " x " << config_ny << " x " << config_nz << ") "
       << "is set.\n";
    std::cout << ss.str();
  }

#endif
}

MpiSupport::MpiGuard::~MpiGuard() {
#if defined(XFDTD_CORE_WITH_MPI)
  auto err = MPI_Finalize();
  if (err != MPI_SUCCESS) {
    std::cerr << "MPI_Finalize failed\n";
    MpiSupport::instance().abort(err);
  }
#endif
}

auto MpiSupport::toString() const -> std::string { return _config.toString(); }

}  // namespace xfdtd
