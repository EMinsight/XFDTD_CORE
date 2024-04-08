#include <xfdtd/common/type_define.h>
#include <xfdtd/parallel/mpi_support.h>
#include <xfdtd/util/fdtd_basic.h>

#include "parallel/mpi_type_define.h"

#if defined(XFDTD_CORE_WITH_MPI)
#include <mpi.h>
#endif

namespace xfdtd {

auto MpiSupport::Slice::makeRowMajorXSlice(std::size_t nx, std::size_t ny,
                                           std::size_t nz)
    -> MpiSupport::Slice {
  Slice slice;
#if defined(XFDTD_CORE_WITH_MPI)
  slice.createVec(1, ny * nz, 1, mpi_type::XFDTD_MPI_REAL_TYPE);
  slice.createStruct(1, {1}, {0});
#endif
  return slice;
}

auto MpiSupport::Slice::makeRowMajorYSlice(std::size_t nx, std::size_t ny,
                                           std::size_t nz)
    -> MpiSupport::Slice {
  Slice slice;
#if defined(XFDTD_CORE_WITH_MPI)
  slice.createVec(1, nz, 1, mpi_type::XFDTD_MPI_REAL_TYPE);
  std::vector<MPI_Aint> offsets(nx, 0);
  for (std::size_t i = 0; i < nx; ++i) {
    offsets[i] = sizeof(Real) * i * ny * nz;
  }
  slice.createStruct(nx, std::vector<int>(nx, 1), std::move(offsets));
#endif
  return slice;
}

auto MpiSupport::Slice::makeRowMajorZSlice(std::size_t nx, std::size_t ny,
                                           std::size_t nz)
    -> MpiSupport::Slice {
  Slice slice;
#if defined(XFDTD_CORE_WITH_MPI)
  slice.createVec(ny, 1, nz, mpi_type::XFDTD_MPI_REAL_TYPE);
  std::vector<MPI_Aint> offsets(nx, 0);
  for (std::size_t i = 0; i < nx; ++i) {
    offsets[i] = sizeof(Real) * i * ny * nz;
  }
  slice.createStruct(nx, std::vector<int>(nx, 1), std::move(offsets));
#endif
  return slice;
}

auto MpiSupport::Slice::makeEyXSlice(std::size_t nx, std::size_t ny,
                                     std::size_t nz) -> MpiSupport::Slice {
  nx = basic::GridStructure::eySizeX(nx);
  ny = basic::GridStructure::eySizeY(ny);
  nz = basic::GridStructure::eySizeZ(nz);
  return makeRowMajorXSlice(nx, ny, nz);
}

auto MpiSupport::Slice::makeEzXSlice(std::size_t nx, std::size_t ny,
                                     std::size_t nz) -> MpiSupport::Slice {
  nx = basic::GridStructure::ezSizeX(nx);
  ny = basic::GridStructure::ezSizeY(ny);
  nz = basic::GridStructure::ezSizeZ(nz);
  return makeRowMajorXSlice(nx, ny, nz);
}

auto MpiSupport::Slice::makeEzYSlice(std::size_t nx, std::size_t ny,
                                     std::size_t nz) -> MpiSupport::Slice {
  nx = basic::GridStructure::ezSizeX(nx);
  ny = basic::GridStructure::ezSizeY(ny);
  nz = basic::GridStructure::ezSizeZ(nz);
  return makeRowMajorYSlice(nx, ny, nz);
}

auto MpiSupport::Slice::makeExYSlice(std::size_t nx, std::size_t ny,
                                     std::size_t nz) -> MpiSupport::Slice {
  nx = basic::GridStructure::exSizeX(nx);
  ny = basic::GridStructure::exSizeY(ny);
  nz = basic::GridStructure::exSizeZ(nz);
  return makeRowMajorYSlice(nx, ny, nz);
}

auto MpiSupport::Slice::makeExZSlice(std::size_t nx, std::size_t ny,
                                     std::size_t nz) -> MpiSupport::Slice {
  nx = basic::GridStructure::exSizeX(nx);
  ny = basic::GridStructure::exSizeY(ny);
  nz = basic::GridStructure::exSizeZ(nz);
  return makeRowMajorZSlice(nx, ny, nz);
}

auto MpiSupport::Slice::makeEyZSlice(std::size_t nx, std::size_t ny,
                                     std::size_t nz) -> MpiSupport::Slice {
  nx = basic::GridStructure::eySizeX(nx);
  ny = basic::GridStructure::eySizeY(ny);
  nz = basic::GridStructure::eySizeZ(nz);
  return makeRowMajorZSlice(nx, ny, nz);
}

auto MpiSupport::Slice::makeHyXSlice(std::size_t nx, std::size_t ny,
                                     std::size_t nz) -> MpiSupport::Slice {
  nx = basic::GridStructure::hySizeX(nx);
  ny = basic::GridStructure::hySizeY(ny);
  nz = basic::GridStructure::hySizeZ(nz);
  return makeRowMajorXSlice(nx, ny, nz);
}

auto MpiSupport::Slice::makeHzXSlice(std::size_t nx, std::size_t ny,
                                     std::size_t nz) -> MpiSupport::Slice {
  nx = basic::GridStructure::hzSizeX(nx);
  ny = basic::GridStructure::hzSizeY(ny);
  nz = basic::GridStructure::hzSizeZ(nz);
  return makeRowMajorXSlice(nx, ny, nz);
}

auto MpiSupport::Slice::makeHzYSlice(std::size_t nx, std::size_t ny,
                                     std::size_t nz) -> MpiSupport::Slice {
  nx = basic::GridStructure::hzSizeX(nx);
  ny = basic::GridStructure::hzSizeY(ny);
  nz = basic::GridStructure::hzSizeZ(nz);
  return makeRowMajorYSlice(nx, ny, nz);
}

auto MpiSupport::Slice::makeHxYSlice(std::size_t nx, std::size_t ny,
                                     std::size_t nz) -> MpiSupport::Slice {
  nx = basic::GridStructure::hxSizeX(nx);
  ny = basic::GridStructure::hxSizeY(ny);
  nz = basic::GridStructure::hxSizeZ(nz);
  return makeRowMajorYSlice(nx, ny, nz);
}

auto MpiSupport::Slice::makeHxZSlice(std::size_t nx, std::size_t ny,
                                     std::size_t nz) -> MpiSupport::Slice {
  nx = basic::GridStructure::hxSizeX(nx);
  ny = basic::GridStructure::hxSizeY(ny);
  nz = basic::GridStructure::hxSizeZ(nz);
  return makeRowMajorZSlice(nx, ny, nz);
}

auto MpiSupport::Slice::makeHyZSlice(std::size_t nx, std::size_t ny,
                                     std::size_t nz) -> MpiSupport::Slice {
  nx = basic::GridStructure::hySizeX(nx);
  ny = basic::GridStructure::hySizeY(ny);
  nz = basic::GridStructure::hySizeZ(nz);
  return makeRowMajorZSlice(nx, ny, nz);
}

#if defined(XFDTD_CORE_WITH_MPI)
auto MpiSupport::Slice::createVec(int count, int block_len, int stripe,
                                  MPI_Datatype type) -> void {
  MPI_Type_vector(count, block_len, stripe, type, &_vec_type);
  MPI_Type_commit(&_vec_type);
}
#endif

#if defined(XFDTD_CORE_WITH_MPI)
auto MpiSupport::Slice::createStruct(int count, std::vector<int> block_lens,
                                     std::vector<MPI_Aint> offsets) -> void {
  _block_lens = std::move(block_lens);
  _offsets = std::move(offsets);
  _vec_types = std::vector<MPI_Datatype>(_block_lens.size(), _vec_type);
  MPI_Type_create_struct(count, _block_lens.data(), _offsets.data(),
                         _vec_types.data(), &_slice);
  MPI_Type_commit(&_slice);
}
#endif

#if defined(XFDTD_CORE_WITH_MPI)
MpiSupport::Slice::Slice(Slice&& other) noexcept
    : _vec_type(other._vec_type),
      _vec_types(std::move(other._vec_types)),
      _block_lens(std::move(other._block_lens)),
      _offsets(std::move(other._offsets)),
      _slice(other._slice) {
  other._vec_type = MPI_DATATYPE_NULL;
  other._slice = MPI_DATATYPE_NULL;
}
#else
MpiSupport::Slice::Slice(Slice&& other) noexcept = default;
#endif

#if defined(XFDTD_CORE_WITH_MPI)
auto MpiSupport::Slice::operator=(Slice&& other) noexcept -> Slice& {
  if (this == &other) {
    return *this;
  }

  if (_vec_type != MPI_DATATYPE_NULL) {
    MPI_Type_free(&_vec_type);
  }

  if (_slice != MPI_DATATYPE_NULL) {
    MPI_Type_free(&_slice);
  }

  _vec_type = other._vec_type;
  _vec_types = std::move(other._vec_types);
  _block_lens = std::move(other._block_lens);
  _offsets = std::move(other._offsets);
  _slice = other._slice;

  other._vec_type = MPI_DATATYPE_NULL;
  other._slice = MPI_DATATYPE_NULL;

  return *this;
}
#else
auto MpiSupport::Slice::operator=(Slice&& other) noexcept -> Slice& = default;
#endif
#if defined(XFDTD_CORE_WITH_MPI)
MpiSupport::Slice::~Slice() {
  if (_vec_type != MPI_DATATYPE_NULL) {
    MPI_Type_free(&_vec_type);
  }

  if (_slice != MPI_DATATYPE_NULL) {
    MPI_Type_free(&_slice);
  }
}
#else
MpiSupport::Slice::~Slice() {}
#endif

}  // namespace xfdtd
