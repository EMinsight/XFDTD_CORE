#ifndef _XFDTD_LIB_UPDATE_SCHEME_H_
#define _XFDTD_LIB_UPDATE_SCHEME_H_

#include <tuple>

namespace xfdtd {

/**
 * FDTD Scheme for updating the electric field.
 * `a,b,c` means the axis.
 * `c` = `a` cross `b`.
 * `p,q` means the different direction.
 * `q` is the behind direction of `p`.
 *
 * @param cece The coefficient cece.
 * @param ec The current value of the electric field.
 * @param cecha The coefficient cecha.
 * @param ha_p The current value of the magnetic field ha_p.
 * @param ha_q The current value of the magnetic field ha_q.
 * @param cechb The coefficient cechb.
 * @param hb_p The current value of the magnetic field hb_p.
 * @param hb_q The current value of the magnetic field hb_q.
 * @return The next time step value of the electric field.
 */
template <typename T>
inline auto eNext(const T& cece, T ec, const T& cecha, T ha_p, T ha_q,
                  const T& cechb, T hb_p, T hb_q) {
  auto&& result = cece * ec + cecha * (ha_p - ha_q) + cechb * (hb_p - hb_q);
  return result;
}

/**
 * FDTD Scheme for updating the magnetic field.
 * `a,b,c` means the axis.
 * `c` = `a` cross `b`.
 * `p,q` means the different direction.
 * `q` is the behind direction of `p`.
 *
 * @param chch The coefficient chch.
 * @param hc The current value of the magnetic field.
 * @param chcea The coefficient chcea.
 * @param ea_p The current value of the electric field ea_p.
 * @param ea_q The current value of the electric field ea_q.
 * @param chceb The coefficient chceb.
 * @param eb_p The current value of the electric field eb_p.
 * @param eb_q The current value of the electric field eb_q.
 * @return The next time step value of the magnetic field.
 */
template <typename T>
inline auto hNext(const T& chch, T hc, const T& chcea, T ea_p, T ea_q,
                  const T& chceb, T eb_p, T eb_q) {
  auto&& result = chch * hc + chcea * (ea_p - ea_q) + chceb * (eb_p - eb_q);
  return result;
}

/**
 * @brief Calculates the shape of a node in a 3D grid.
 *
 * This function takes the start and end indices of each dimension (i, j, k) and
 * calculates the number of nodes in each dimension (nx, ny, nz).
 *
 * @param i_start The start index of dimension i.
 * @param i_end The end index of dimension i.
 * @param j_start The start index of dimension j.
 * @param j_end The end index of dimension j.
 * @param k_start The start index of dimension k.
 * @param k_end The end index of dimension k.
 * @return A tuple containing the number of nodes in each dimension (nx, ny,
 * nz).
 */
template <typename T>
inline auto exNodeShape(const T& i_start, const T& i_end, const T& j_start,
                        const T& j_end, const T& k_start, const T& k_end) {
  auto&& nx = i_end - i_start;
  auto&& ny = j_end - j_start + 1;
  auto&& nz = k_end - k_start + 1;
  return std::make_tuple(nx, ny, nz);
}

/**
 * @brief Calculates the shape of a node in a 3D grid.
 *
 * This function takes the start and end indices of each dimension (i, j, k) and
 * calculates the number of nodes in each dimension (nx, ny, nz).
 *
 * @param i_start The start index of dimension i.
 * @param i_end The end index of dimension i.
 * @param j_start The start index of dimension j.
 * @param j_end The end index of dimension j.
 * @param k_start The start index of dimension k.
 * @param k_end The end index of dimension k.
 * @return A tuple containing the number of nodes in each dimension (nx, ny,
 * nz).
 */
template <typename T>
inline auto eyNodeShape(const T& i_start, const T& i_end, const T& j_start,
                        const T& j_end, const T& k_start, const T& k_end) {
  auto&& nx = i_end - i_start + 1;
  auto&& ny = j_end - j_start;
  auto&& nz = k_end - k_start + 1;
  return std::make_tuple(nx, ny, nz);
}

/**
 * @brief Calculates the shape of a node in a 3D grid.
 *
 * This function takes the start and end indices of each dimension (i, j, k) and
 * calculates the number of nodes in each dimension (nx, ny, nz).
 *
 * @param i_start The start index of dimension i.
 * @param i_end The end index of dimension i.
 * @param j_start The start index of dimension j.
 * @param j_end The end index of dimension j.
 * @param k_start The start index of dimension k.
 * @param k_end The end index of dimension k.
 * @return A tuple containing the number of nodes in each dimension (nx, ny,
 * nz).
 */
template <typename T>
inline auto ezNodeShape(const T& i_start, const T& i_end, const T& j_start,
                        const T& j_end, const T& k_start, const T& k_end) {
  auto&& nx = i_end - i_start + 1;
  auto&& ny = j_end - j_start + 1;
  auto&& nz = k_end - k_start;
  return std::make_tuple(nx, ny, nz);
}

/**
 * @brief Calculates the shape of a node in a 3D grid.
 *
 * This function takes the start and end indices of each dimension (i, j, k) and
 * calculates the number of nodes in each dimension (nx, ny, nz).
 *
 * @param i_start The start index of dimension i.
 * @param i_end The end index of dimension i.
 * @param j_start The start index of dimension j.
 * @param j_end The end index of dimension j.
 * @param k_start The start index of dimension k.
 * @param k_end The end index of dimension k.
 * @return A tuple containing the number of nodes in each dimension (nx, ny,
 * nz).
 */
template <typename T>
inline auto hxNodeShape(const T& i_start, const T& i_end, const T& j_start,
                        const T& j_end, const T& k_start, const T& k_end) {
  auto&& nx = i_end - i_start + 1;
  auto&& ny = j_end - j_start;
  auto&& nz = k_end - k_start;
  return std::make_tuple(nx, ny, nz);
}

/**
 * @brief Calculates the shape of a node in a 3D grid.
 *
 * This function takes the start and end indices of each dimension (i, j, k) and
 * calculates the number of nodes in each dimension (nx, ny, nz).
 *
 * @param i_start The start index of dimension i.
 * @param i_end The end index of dimension i.
 * @param j_start The start index of dimension j.
 * @param j_end The end index of dimension j.
 * @param k_start The start index of dimension k.
 * @param k_end The end index of dimension k.
 * @return A tuple containing the number of nodes in each dimension (nx, ny,
 * nz).
 */
template <typename T>
inline auto hyNodeShape(const T& i_start, const T& i_end, const T& j_start,
                        const T& j_end, const T& k_start, const T& k_end) {
  auto&& nx = i_end - i_start;
  auto&& ny = j_end - j_start + 1;
  auto&& nz = k_end - k_start;
  return std::make_tuple(nx, ny, nz);
}

/**
 * @brief Calculates the shape of a node in a 3D grid.
 *
 * This function takes the start and end indices of each dimension (i, j, k) and
 * calculates the number of nodes in each dimension (nx, ny, nz).
 *
 * @param i_start The start index of dimension i.
 * @param i_end The end index of dimension i.
 * @param j_start The start index of dimension j.
 * @param j_end The end index of dimension j.
 * @param k_start The start index of dimension k.
 * @param k_end The end index of dimension k.
 * @return A tuple containing the number of nodes in each dimension (nx, ny,
 * nz).
 */
template <typename T>
inline auto hzNodeShape(const T& i_start, const T& i_end, const T& j_start,
                        const T& j_end, const T& k_start, const T& k_end) {
  auto&& nx = i_end - i_start;
  auto&& ny = j_end - j_start;
  auto&& nz = k_end - k_start + 1;
  return std::make_tuple(nx, ny, nz);
}

template <typename Index, typename T>
inline auto updateEzEdgeYZ(const Index i, const Index js, const Index je,
                           const Index ks, const Index ke, const T& ceze,
                           const T& cezhx, const T& cezhy, const T& hx,
                           const T& hy, const T& hy_buffer, T& ez) {
  for (Index j{js + 1}; j < je; ++j) {
    for (Index k{ks}; k < ke; ++k) {
      ez(i, j, k) = eNext(ceze(i, j, k), ez(i, j, k), cezhx(i, j, k),
                          hx(i, j, k), hx(i, j - 1, k), cezhy(i, j, k),
                          hy(i, j, k), hy_buffer(0, (j - js), k - ks));
    }
  }
}

template <typename Index, typename T>
inline auto updateEzEdgeXZ(const Index is, const Index ie, const Index j,
                           const Index ks, const Index ke, const T& ceze,
                           const T& cezhx, const T& cezhy, const T& hx,
                           const T& hx_buffer, const T& hy, T& ez) {
  for (std::size_t i{is + 1}; i < ie; ++i) {
    for (std::size_t k{ks}; k < ke; ++k) {
      ez(i, j, k) = eNext(ceze(i, j, k), ez(i, j, k), cezhx(i, j, k),
                          hx(i, j, k), hx_buffer(i - is, 0, k - ks),
                          cezhy(i, j, k), hy(i, j, k), hy(i - 1, j, k));
    }
  }
}

template <typename Index, typename T>
inline auto updateEzCornerXY(const Index i, const Index j, const Index ks,
                             const Index ke, const T& ceze, const T& cezhx,
                             const T& cezhy, const T& hx, const T& hx_buffer,
                             const T& hy, const T& hy_buffer, T& ez) {
  for (std::size_t k{ks}; k < ke; ++k) {
    ez(i, j, k) = eNext(ceze(i, j, k), ez(i, j, k), cezhx(i, j, k), hx(i, j, k),
                        hx_buffer(0, 0, k - ks), cezhy(i, j, k), hy(i, j, k),
                        hy_buffer(0, 0, k - ks));
  }
}

}  // namespace xfdtd

#endif  // _XFDTD_LIB_UPDATE_SCHEME_H_
