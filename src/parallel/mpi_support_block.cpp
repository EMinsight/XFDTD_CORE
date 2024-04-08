#include <xfdtd/common/type_define.h>
#include <xfdtd/parallel/mpi_support.h>

#include "parallel/mpi_type_define.h"

namespace xfdtd {

auto MpiSupport::Block::make(Profile profile) -> Block {
  Block block;

  block._profile = profile;

#if defined(XFDTD_CORE_WITH_MPI)
  auto nx = profile._nx;
  auto ny = profile._ny;
  auto nz = profile._nz;
  auto stride_vec = profile._stride_vec;
  auto stride_elem = profile._stride_elem;
  MPI_Type_vector(ny, nz, stride_elem, mpi_type::XFDTD_MPI_REAL_TYPE,
                  &block._vec_type);
  MPI_Type_commit(&block._vec_type);
  block._vec_types = std::vector<MPI_Datatype>(nx, block._vec_type);
  block._block_lens = std::vector<int>(nx, 1);
  block._offsets = std::vector<MPI_Aint>(nx, 0);
  for (int i = 0; i < nx; ++i) {
    block._offsets[i] = sizeof(Real) * i * stride_vec;
  }
  MPI_Type_create_struct(nx, block._block_lens.data(), block._offsets.data(),
                         block._vec_types.data(), &block._block);
  MPI_Type_commit(&block._block);
#endif

  return block;
}

#if defined(XFDTD_CORE_WITH_MPI)

MpiSupport::Block::Block(Block&& other) noexcept
    : _profile{other._profile},
      _vec_type{other._vec_type},
      _vec_types{std::move(other._vec_types)},
      _block_lens{std::move(other._block_lens)},
      _offsets{std::move(other._offsets)},
      _block{other._block} {
  other._profile = Block::Profile{};
  other._vec_type = MPI_DATATYPE_NULL;
  other._block = MPI_DATATYPE_NULL;
}

auto MpiSupport::Block::operator=(Block&& other) noexcept -> Block& {
  if (this != &other) {
    _profile = other._profile;
    _vec_type = other._vec_type;
    _vec_types = std::move(other._vec_types);
    _block_lens = std::move(other._block_lens);
    _offsets = std::move(other._offsets);
    _block = other._block;

    other._profile = Block::Profile{};
    other._vec_type = MPI_DATATYPE_NULL;
    other._block = MPI_DATATYPE_NULL;
  }

  return *this;
}

MpiSupport::Block::~Block() {
  if (_vec_type != MPI_DATATYPE_NULL) {
    MPI_Type_free(&_vec_type);
  }
  if (_block != MPI_DATATYPE_NULL) {
    MPI_Type_free(&_block);
  }
}

#else

MpiSupport::Block::Block(Block&& other) noexcept : _profile{other._profile} {}

auto MpiSupport::Block::operator=(Block&& other) noexcept -> Block& {
  if (this != &other) {
    _profile = other._profile;
  }

  return *this;
}

MpiSupport::Block::~Block() {}

#endif

auto MpiSupport::Block::profile() const -> const Profile& { return _profile; }

}  // namespace xfdtd
