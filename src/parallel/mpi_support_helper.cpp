/**
 * @file mpi_support_helper.cpp
 * @author your name (you@domain.com)
 * @brief Wrapper for MPI functions
 * @version 0.1
 * @date 2024-03-27
 *
 * @copyright Copyright (c) 2024
 *
 */
#include <xfdtd/parallel/mpi_support.h>

#include "parallel/mpi_type_define.h"

#if defined(XFDTD_CORE_WITH_MPI)
#include <mpi.h>
#endif

namespace xfdtd {

auto MpiSupport::send(const MpiConfig& config, const void* buf, int count,
                      int dest, int tag) -> void {
#if defined(XFDTD_CORE_WITH_MPI)
  MPI_Send(buf, count, MPI_BYTE, dest, tag, config.comm());
#endif
}

auto MpiSupport::send(const MpiConfig& config, const void* buf, int count,
                      const TypeGuard& type, int dest, int tag) -> void {
#if defined(XFDTD_CORE_WITH_MPI)
  MPI_Send(buf, count, type._type, dest, tag, config.comm());
#endif
}

auto MpiSupport::send(const MpiConfig& config, const Real* buf, int count,
                      const Block& block, int dest, int tag) -> void {
#if defined(XFDTD_CORE_WITH_MPI)
  MPI_Send(&buf[block.profile()._disp], count, block.block(), dest, tag,
           config.comm());
#endif
}

auto MpiSupport::iSend(const MpiConfig& config, const void* buf, int count,
                       int dest, int tag) -> std::size_t {
#if defined(XFDTD_CORE_WITH_MPI)
  MPI_Request request;
  MPI_Isend(buf, count, MPI_BYTE, dest, tag, config.comm(), &request);
  _requests.emplace_back(request);
  return _requests.size() - 1;
#endif
  return -1;
}

auto MpiSupport::iSend(const MpiConfig& config, const void* buf, int count,
                       const TypeGuard& type, int dest, int tag)
    -> std::size_t {
#if defined(XFDTD_CORE_WITH_MPI)
  MPI_Request request;
  MPI_Isend(buf, count, type._type, dest, tag, config.comm(), &request);
  _requests.emplace_back(request);
  return _requests.size() - 1;
#endif
  return -1;
}

auto MpiSupport::iSend(const MpiConfig& config, const Real* buf, int count,
                       const Block& block, int dest, int tag) -> std::size_t {
#if defined(XFDTD_CORE_WITH_MPI)
  MPI_Request request;
  MPI_Isend(&buf[block.profile()._disp], count, block.block(), dest, tag,
            config.comm(), &request);
  _requests.emplace_back(request);
  return _requests.size() - 1;
#endif
  return -1;
}

auto MpiSupport::recv(const MpiConfig& config, void* buf, int count, int source,
                      int tag) -> void {
#if defined(XFDTD_CORE_WITH_MPI)
  MPI_Recv(buf, count, MPI_BYTE, source, tag, config.comm(), MPI_STATUS_IGNORE);
#endif
}

auto MpiSupport::recv(const MpiConfig& config, void* buf, int count,
                      const TypeGuard& type, int source, int tag) -> void {
#if defined(XFDTD_CORE_WITH_MPI)
  MPI_Recv(buf, count, type._type, source, tag, config.comm(),
           MPI_STATUS_IGNORE);
#endif
}

auto MpiSupport::recv(const MpiConfig& config, Real* buf, int count,
                      const Block& block, int source, int tag) -> void {
#if defined(XFDTD_CORE_WITH_MPI)
  MPI_Recv(&buf[block.profile()._disp], count, block.block(), source, tag,
           config.comm(), MPI_STATUS_IGNORE);
#endif
}

auto MpiSupport::iRecv(const MpiConfig& config, void* buf, int count,
                       int source, int tag) -> std::size_t {
#if defined(XFDTD_CORE_WITH_MPI)
  MPI_Request request;
  MPI_Irecv(buf, count, MPI_BYTE, source, tag, config.comm(), &request);
  _requests.emplace_back(request);
  return _requests.size() - 1;
#endif
  return -1;
}

auto MpiSupport::iRecv(const MpiConfig& config, void* buf, int count,
                       const TypeGuard& type, int source, int tag)
    -> std::size_t {
#if defined(XFDTD_CORE_WITH_MPI)
  MPI_Request request;
  MPI_Irecv(buf, count, type._type, source, tag, config.comm(), &request);
  _requests.emplace_back(request);
  return _requests.size() - 1;
#endif
  return -1;
}

auto MpiSupport::iRecv(const MpiConfig& config, Real* buf, int count,
                       const Block& block, int source, int tag) -> std::size_t {
#if defined(XFDTD_CORE_WITH_MPI)
  MPI_Request request;
  MPI_Irecv(&buf[block.profile()._disp], count, block.block(), source, tag,
            config.comm(), &request);
  _requests.emplace_back(request);
  return _requests.size() - 1;
#endif
  return -1;
}

auto MpiSupport::sendRecv(const MpiConfig& config, const void* send_buf,
                          int send_count, int dest, int send_tag,
                          void* recv_buf, int recv_count, int source,
                          int recv_tag) -> void {
#if defined(XFDTD_CORE_WITH_MPI)
  MPI_Sendrecv(send_buf, send_count, MPI_BYTE, dest, send_tag, recv_buf,
               recv_count, MPI_BYTE, source, recv_tag, config.comm(),
               MPI_STATUS_IGNORE);
#endif
}

auto MpiSupport::sendRecv(const MpiConfig& config, const void* send_buf,
                          int send_count, const TypeGuard& send_type, int dest,
                          int send_tag, void* recv_buf, int recv_count,
                          const TypeGuard& recv_type, int source, int recv_tag)
    -> void {
#if defined(XFDTD_CORE_WITH_MPI)
  MPI_Sendrecv(send_buf, send_count, send_type._type, dest, send_tag, recv_buf,
               recv_count, recv_type._type, source, recv_tag, config.comm(),
               MPI_STATUS_IGNORE);
#endif
}

auto MpiSupport::sendRecv(const MpiConfig& config, const Real* send_buf,
                          int send_count, int dest, int send_tag,
                          Real* recv_buf, int recv_count,
                          const Block& recv_block, int source, int recv_tag)
    -> void {
#if defined(XFDTD_CORE_WITH_MPI)
  MPI_Sendrecv(send_buf, send_count, mpi_type::XFDTD_MPI_REAL_TYPE, dest,
               send_tag, &recv_buf[recv_block.profile()._disp], recv_count,
               recv_block.block(), source, recv_tag, config.comm(),
               MPI_STATUS_IGNORE);
#endif
}

auto MpiSupport::iSendRecv(const MpiConfig& config, const void* send_buf,
                           int send_count, int dest, int send_tag,
                           void* recv_buf, int recv_count, int source,
                           int recv_tag) -> std::size_t {
#if defined(XFDTD_CORE_WITH_MPI)
  MPI_Request request;
  MPI_Isendrecv(send_buf, send_count, MPI_BYTE, dest, send_tag, recv_buf,
                recv_count, MPI_BYTE, source, recv_tag, config.comm(),
                &request);
  _requests.emplace_back(request);
  return _requests.size() - 1;
#endif
  return -1;
}

auto MpiSupport::iSendRecv(const MpiConfig& config, const void* send_buf,
                           int send_count, const TypeGuard& send_type, int dest,
                           int send_tag, void* recv_buf, int recv_count,
                           const TypeGuard& recv_type, int source, int recv_tag)
    -> std::size_t {
#if defined(XFDTD_CORE_WITH_MPI)
  MPI_Request request;
  MPI_Isendrecv(send_buf, send_count, send_type._type, dest, send_tag, recv_buf,
                recv_count, recv_type._type, source, recv_tag, config.comm(),
                &request);
  _requests.emplace_back(request);
  return _requests.size() - 1;
#endif
  return -1;
}

auto MpiSupport::iSendRecv(const MpiConfig& config, const Real* send_buf,
                           int send_count, int dest, int send_tag,
                           Real* recv_buf, int recv_count,
                           const Block& recv_block, int source, int recv_tag)
    -> std::size_t {
#if defined(XFDTD_CORE_WITH_MPI)
  MPI_Request request;
  MPI_Isendrecv(send_buf, send_count, mpi_type::XFDTD_MPI_REAL_TYPE, dest,
                send_tag, &recv_buf[recv_block.profile()._disp], recv_count,
                recv_block.block(), source, recv_tag, config.comm(), &request);
  _requests.emplace_back(request);
  return _requests.size() - 1;
#endif
  return -1;
}

auto MpiSupport::gather(const MpiConfig& config, const void* send_buf,
                        int send_count, void* recv_buf, int recv_count,
                        int root) -> void {
#if defined(XFDTD_CORE_WITH_MPI)
  MPI_Gather(send_buf, send_count, MPI_BYTE, recv_buf, recv_count, MPI_BYTE,
             root, config.comm());
#endif
}

auto MpiSupport::gather(const MpiConfig& config, const void* send_buf,
                        int send_count, const TypeGuard& send_type,
                        void* recv_buf, int recv_count,
                        const TypeGuard& recv_type, int root) -> void {
#if defined(XFDTD_CORE_WITH_MPI)
  MPI_Gather(send_buf, send_count, send_type._type, recv_buf, recv_count,
             recv_type._type, root, config.comm());
#endif
}

auto MpiSupport::gather(const MpiConfig& config, const Real* send_buf,
                        int send_count, Real* recv_buf, int recv_count,
                        const Block& recv_block, int root) -> void {
#if defined(XFDTD_CORE_WITH_MPI)
  MPI_Gather(send_buf, send_count, recv_block.block(),
             &recv_buf[recv_block.profile()._disp], recv_count,
             recv_block.block(), root, config.comm());
#endif
}

auto MpiSupport::iGather(const MpiConfig& config, const void* send_buf,
                         int send_count, void* recv_buf, int recv_count,
                         int root) -> std::size_t {
#if defined(XFDTD_CORE_WITH_MPI)
  MPI_Request request;
  MPI_Igather(send_buf, send_count, MPI_BYTE, recv_buf, recv_count, MPI_BYTE,
              root, config.comm(), &request);
  _requests.emplace_back(request);
  return _requests.size() - 1;
#endif
  return -1;
}

auto MpiSupport::iGather(const MpiConfig& config, const void* send_buf,
                         int send_count, const TypeGuard& send_type,
                         void* recv_buf, int recv_count,
                         const TypeGuard& recv_type, int root) -> std::size_t {
#if defined(XFDTD_CORE_WITH_MPI)
  MPI_Request request;
  MPI_Igather(send_buf, send_count, send_type._type, recv_buf, recv_count,
              recv_type._type, root, config.comm(), &request);
  _requests.emplace_back(request);
  return _requests.size() - 1;
#endif
  return -1;
}

auto MpiSupport::iGather(const MpiConfig& config, const Real* send_buf,
                         int send_count, Real* recv_buf, int recv_count,
                         const Block& recv_block, int root) -> std::size_t {
#if defined(XFDTD_CORE_WITH_MPI)
  MPI_Request request;
  MPI_Igather(send_buf, send_count, recv_block.block(),
              &recv_buf[recv_block.profile()._disp], recv_count,
              recv_block.block(), root, config.comm(), &request);
  _requests.emplace_back(request);
  return _requests.size() - 1;
#endif
  return -1;
}

auto MpiSupport::allGather(const MpiConfig& config, const void* send_buf,
                           int send_count, void* recv_buf, int recv_count)
    -> void {
#if defined(XFDTD_CORE_WITH_MPI)
  MPI_Allgather(send_buf, send_count, MPI_BYTE, recv_buf, recv_count, MPI_BYTE,
                config.comm());
#endif
}

auto MpiSupport::reduceSum(const MpiConfig& config, const Real* send_buf,
                           Real* recv_buf, int count) -> void {
#if defined(XFDTD_CORE_WITH_MPI)
  MPI_Reduce(send_buf, recv_buf, count, mpi_type::XFDTD_MPI_REAL_TYPE, MPI_SUM,
             config.root(), config.comm());
#endif
}

auto MpiSupport::reduceSum(const MpiConfig& config,
                           const std::complex<Real>* send_buf,
                           std::complex<Real>* recv_buf, int count) const
    -> void {
#if defined(XFDTD_CORE_WITH_MPI)
#if defined(XFDTD_CORE_SINGLE_PRECISION)
  MPI_Reduce(send_buf, recv_buf, count, MPI_COMPLEX, MPI_SUM, config.root(),
             config.comm());
#else
  MPI_Reduce(send_buf, recv_buf, count, mpi_type::XFDTD_MPI_COMPLEX_TYPE,
             MPI_SUM, config.root(), config.comm());
#endif
#endif
}

}  // namespace xfdtd
