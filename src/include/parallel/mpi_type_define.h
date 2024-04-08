#ifndef __XFDTD_CORE_MPI_TYPE_DEFINE_H__
#define __XFDTD_CORE_MPI_TYPE_DEFINE_H__

#if defined(XFDTD_CORE_WITH_MPI)
#include <mpi.h>
#endif

namespace xfdtd::mpi_type {
#if defined(XFDTD_CORE_WITH_MPI)

#if defined(XFDTD_CORE_SINGLE_PRECISION)
static const auto XFDTD_MPI_REAL_TYPE = MPI_FLOAT;
static const auto XFDTD_MPI_COMPLEX_TYPE = MPI_COMPLEX;
#else
static const auto XFDTD_MPI_REAL_TYPE = MPI_DOUBLE;
static const auto XFDTD_MPI_COMPLEX_TYPE = MPI_DOUBLE_COMPLEX;
#endif

#endif
}  // namespace xfdtd::mpi_type

#endif  // __XFDTD_CORE_MPI_TYPE_DEFINE_H__
