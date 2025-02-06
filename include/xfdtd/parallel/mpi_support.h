#ifndef __XFDTD_CORE_MPI_SUPPORT_H__
#define __XFDTD_CORE_MPI_SUPPORT_H__

#include <xfdtd/common/type_define.h>
#include <xfdtd/exception/exception.h>
#include <xfdtd/parallel/mpi_config.h>

#include <complex>

#if defined(XFDTD_CORE_WITH_MPI)
#include <mpi.h>
#endif

namespace xfdtd {

class XFDTDMpiSupportException : public XFDTDException {
 public:
  explicit XFDTDMpiSupportException(const std::string& message = "")
      : XFDTDException{message} {}
};

class MpiSupport {
 public:
  /**
   * @brief prepare MPI environment
   *
   * @return true If MPI is initialized successfully or already initialized.
   * @return false If MPI is failed to initialize or not compiled with MPI.
   */
  static auto init(int argc = 0, char** argv = nullptr) -> bool;

  /**
   * @brief Set the dimension of the MPI parallel layout. It will change the
   * configuration of MPI whether it returns true or false. It's not valid to
   * call this after geting the MpiSupport instance.
   *
   * @return true if the dimension is matched with the number of processes.
   */
  static auto setMpiParallelDim(int nx, int ny, int nz) -> bool;

  static auto globalRank() -> int;

  static auto globalSize() -> int;

 public:
  class Block {
   public:
    struct Profile {
      int _nx{1};
      int _ny{1};
      int _nz{0};
      int _stride_vec{1};
      int _stride_elem{1};
      int _disp{0};
    };

   public:
    static auto makeRowMajorXSlice(std::size_t thickness_x, std::size_t ny,
                                   std::size_t nz, std::size_t disp) -> Block;

    static auto makeRowMajorYSlice(std::size_t thickness_y, std::size_t nx,
                                   std::size_t ny, std::size_t nz,
                                   std::size_t disp) -> Block;

    static auto makeRowMajorZSlice(std::size_t thickness_z, std::size_t nx,
                                   std::size_t ny, std::size_t nz,
                                   std::size_t disp) -> Block;

    static auto make(Profile profile) -> Block;

   public:
    Block() = default;

    Block(Block const&) = delete;

    Block(Block&& other) noexcept;

    auto operator=(Block const&) -> Block& = delete;

    auto operator=(Block&& other) noexcept -> Block&;

    ~Block();

    auto profile() const -> const Profile&;

#if defined(XFDTD_CORE_WITH_MPI)
    auto block() const -> MPI_Datatype { return _block; }
#endif

   private:
    Profile _profile;

#if defined(XFDTD_CORE_WITH_MPI)
    MPI_Datatype _vec_type{MPI_DATATYPE_NULL};

    std::vector<MPI_Datatype> _vec_types;
    std::vector<int> _block_lens;
    std::vector<MPI_Aint> _offsets;

    MPI_Datatype _block{MPI_DATATYPE_NULL};
#endif
  };

  struct TypeGuard {
    explicit TypeGuard() = default;
    TypeGuard(const TypeGuard&) = delete;
    TypeGuard(TypeGuard&&) = delete;
    TypeGuard& operator=(const TypeGuard&) = delete;
    TypeGuard& operator=(TypeGuard&&) = delete;
    ~TypeGuard() {
#if defined(XFDTD_CORE_WITH_MPI)
      if (_type != MPI_DATATYPE_NULL) {
        auto err = MPI_Type_free(&_type);
        if (err != MPI_SUCCESS) {
          std::cerr
              << "MpiSupport::TypeGuard::~TypeGuard MPI_Type_free failed\n";
          MpiSupport::instance().abort(err);
        }
      }
#endif
    }
#if defined(XFDTD_CORE_WITH_MPI)
    MPI_Datatype _type{MPI_DATATYPE_NULL};
#else
    void* _type{nullptr};
#endif
  };

 public:
#if defined(XFDTD_CORE_WITH_MPI)
  inline static constexpr int ANY_SOURCE = MPI_ANY_SOURCE;
  inline static constexpr int ANY_TAG = MPI_ANY_TAG;
#else
  inline static constexpr int ANY_SOURCE = -1;
  inline static constexpr int ANY_TAG = -1;
#endif

  inline static constexpr int EXCHANGE_HY_X_SR_TAG = 0;
  inline static constexpr int EXCHANGE_HY_X_RS_TAG = 1;
  inline static constexpr int EXCHANGE_HZ_X_SR_TAG = 2;
  inline static constexpr int EXCHANGE_HZ_X_RS_TAG = 3;

  inline static constexpr int EXCHANGE_HZ_Y_SR_TAG = 4;
  inline static constexpr int EXCHANGE_HZ_Y_RS_TAG = 5;
  inline static constexpr int EXCHANGE_HX_Y_SR_TAG = 6;
  inline static constexpr int EXCHANGE_HX_Y_RS_TAG = 7;

  inline static constexpr int EXCHANGE_HX_Z_SR_TAG = 8;
  inline static constexpr int EXCHANGE_HX_Z_RS_TAG = 9;
  inline static constexpr int EXCHANGE_HY_Z_SR_TAG = 10;
  inline static constexpr int EXCHANGE_HY_Z_RS_TAG = 11;

  static auto instance(int argc = 0, char** argv = nullptr) -> MpiSupport&;

 public:
  MpiSupport(MpiSupport const&) = delete;

  MpiSupport(MpiSupport&&) = delete;

  auto operator=(MpiSupport const&) -> MpiSupport& = delete;

  auto operator=(MpiSupport&&) -> MpiSupport& = delete;

  ~MpiSupport();

  auto serializedPrint(std::string_view message) -> void;

  auto rank() const -> int { return _config.rank(); }

  auto size() const -> int { return _config.size(); }

  auto root() const -> int { return _config.root(); }

  auto isRoot() const -> bool { return _config.root() == _config.rank(); }

  auto abort(int error_code) const -> void;

  auto barrier() const -> void;

  auto generateSlice(std::size_t nx, std::size_t ny, std::size_t nz) -> void;

  auto config() const -> const MpiConfig& { return _config; }

  auto xNext() const -> int;

  auto xPrev() const -> int;

  auto yNext() const -> int;

  auto yPrev() const -> int;

  auto zNext() const -> int;

  auto zPrev() const -> int;

  auto sendRecvHyXHead(Array3D<Real>& hy) -> void;

  auto recvSendHyXTail(Array3D<Real>& hy) -> void;

  auto sendRecvHzXHead(Array3D<Real>& hz) -> void;

  auto recvSendHzXTail(Array3D<Real>& hz) -> void;

  auto sendRecvHzYHead(Array3D<Real>& hz) -> void;

  auto recvSendHzYTail(Array3D<Real>& hz) -> void;

  auto sendRecvHxYHead(Array3D<Real>& hx) -> void;

  auto recvSendHxYTail(Array3D<Real>& hx) -> void;

  auto sendRecvHxZHead(Array3D<Real>& hx) -> void;

  auto recvSendHxZTail(Array3D<Real>& hx) -> void;

  auto sendRecvHyZHead(Array3D<Real>& hy) -> void;

  auto recvSendHyZTail(Array3D<Real>& hy) -> void;

  auto waitAll() -> void;

  auto toString() const -> std::string;

  /**
   * @brief Wrapper for MPI_Send
   *
   * @param config
   * @param buf
   * @param count unit: byte
   * @param dest
   * @param tag
   */
  auto send(const MpiConfig& config, const void* buf, int count, int dest,
            int tag) -> void;

  auto send(const MpiConfig& config, const void* buf, int count,
            const TypeGuard& type, int dest, int tag) -> void;

  auto send(const MpiConfig& config, const Real* buf, int count,
            const Block& block, int dest, int tag) -> void;

  /**
   * @brief Wrapper for MPI_Isend. add request to the request list. call waitAll
   * to wait for all requests to be completed. Or call wait(request_id) to wait.
   * Unsafe thread.
   *
   * @param config
   * @param buf
   * @param count
   * @param dest
   * @param tag
   * @return std::size_t: request id
   */
  auto iSend(const MpiConfig& config, const void* buf, int count, int dest,
             int tag) -> std::size_t;

  auto iSend(const MpiConfig& config, const void* buf, int count,
             const TypeGuard& type, int dest, int tag) -> std::size_t;

  auto iSend(const MpiConfig& config, const Real* buf, int count,
             const Block& block, int dest, int tag) -> std::size_t;

  auto recv(const MpiConfig& config, void* buf, int count, int source,
            int tag) -> void;

  auto recv(const MpiConfig& config, void* buf, int count,
            const TypeGuard& type, int source, int tag) -> void;

  auto recv(const MpiConfig& config, Real* buf, int count, const Block& block,
            int source, int tag) -> void;

  auto iRecv(const MpiConfig& config, void* buf, int count, int source,
             int tag) -> std::size_t;

  auto iRecv(const MpiConfig& config, void* buf, int count,
             const TypeGuard& type, int source, int tag) -> std::size_t;

  auto iRecv(const MpiConfig& config, Real* buf, int count, const Block& block,
             int source, int tag) -> std::size_t;

  auto sendRecv(const MpiConfig& config, const void* send_buf, int send_count,
                int dest, int send_tag, void* recv_buf, int recv_count,
                int source, int recv_tag) -> void;

  auto sendRecv(const MpiConfig& config, const void* send_buf, int send_count,
                const TypeGuard& send_type, int dest, int send_tag,
                void* recv_buf, int recv_count, const TypeGuard& recv_type,
                int source, int recv_tag) -> void;

  auto sendRecv(const MpiConfig& config, const Real* send_buf, int send_count,
                int dest, int send_tag, Real* recv_buf, int recv_count,
                const Block& recv_block, int source, int recv_tag) -> void;

  auto gather(const MpiConfig& config, const void* send_buf, int send_count,
              void* recv_buf, int recv_count, int root) -> void;

  auto gather(const MpiConfig& config, const void* send_buf, int send_count,
              const TypeGuard& send_type, void* recv_buf, int recv_count,
              const TypeGuard& recv_type, int root) -> void;

  auto gather(const MpiConfig& config, const Real* send_buf, int send_count,
              Real* recv_buf, int recv_count, const Block& recv_block,
              int root) -> void;

  auto iGather(const MpiConfig& config, const void* send_buf, int send_count,
               void* recv_buf, int recv_count, int root) -> std::size_t;

  auto iGather(const MpiConfig& config, const void* send_buf, int send_count,
               const TypeGuard& send_type, void* recv_buf, int recv_count,
               const TypeGuard& recv_type, int root) -> std::size_t;

  auto iGather(const MpiConfig& config, const Real* send_buf, int send_count,
               Real* recv_buf, int recv_count, const Block& recv_block,
               int root) -> std::size_t;

  auto allGather(const MpiConfig& config, const void* send_buf, int send_count,
                 void* recv_buf, int recv_count) -> void;

  auto reduceSum(const MpiConfig& config, const Real* send_buf, Real* recv_buf,
                 int count) -> void;

  auto reduceSum(const MpiConfig& config, const std::complex<Real>* send_buf,
                 std::complex<Real>* recv_buf, int count) const -> void;

 private:
  inline static int global_rank{0};
  inline static int global_size{1};
  inline static int config_nx{1};
  inline static int config_ny{1};
  inline static int config_nz{1};

 private:
  explicit MpiSupport(int argc = 0, char** argv = nullptr);

  struct MpiGuard {
    ~MpiGuard();
  };
  // this is a guard to finalize MPI when the program ends
  MpiGuard _guard{};

  MpiConfig _config;

  // Block must be destroyed before MPI_Finalize.
  Block _hy_xn_r_block;
  Block _hy_xn_s_block;
  Block _hz_xn_r_block;
  Block _hz_xn_s_block;
  Block _hx_yn_r_block;
  Block _hx_yn_s_block;
  Block _hz_yn_r_block;
  Block _hz_yn_s_block;
  Block _hx_zn_r_block;
  Block _hx_zn_s_block;
  Block _hy_zn_r_block;
  Block _hy_zn_s_block;

  Block _hy_xp_r_block;
  Block _hy_xp_s_block;
  Block _hz_xp_r_block;
  Block _hz_xp_s_block;
  Block _hx_yp_r_block;
  Block _hx_yp_s_block;
  Block _hz_yp_r_block;
  Block _hz_yp_s_block;
  Block _hx_zp_r_block;
  Block _hx_zp_s_block;
  Block _hy_zp_r_block;
  Block _hy_zp_s_block;

#if defined(XFDTD_CORE_WITH_MPI)
  std::vector<MPI_Request> _requests;
#endif
};

}  // namespace xfdtd

#endif  // __XFDTD_CORE_MPI_SUPPORT_H__
