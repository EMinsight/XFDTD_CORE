#ifndef __XFDTD_CORE_DECOMPOSE_TASK_H__
#define __XFDTD_CORE_DECOMPOSE_TASK_H__

#include <xfdtd/common/index_task.h>
#include <xfdtd/exception/exception.h>

namespace xfdtd {

template <typename T>
inline auto decomposeSize(int my_rank, int num_procs, const T& problem_size)
    -> Range<T> {
  if (num_procs <= 0) {
    throw XFDTDException("Number of processes is less than or equal to zero");
  }

  if (problem_size < num_procs) {
    throw XFDTDException("Problem size is less than number of processes");
  }

  auto quotient = problem_size / num_procs;
  auto remainder = problem_size % num_procs;
  auto start = T{};
  auto end = T{};

  if (my_rank < remainder) {
    start = my_rank * (quotient + 1);
    end = start + quotient + 1;
    return {start, end};
  }

  start = my_rank * quotient + remainder;
  end = start + quotient;
  return {start, end};
}

template <typename T>
inline auto decomposeRange(const Range<T>& problem, int num_procs)
    -> std::vector<Range<T>> {
  std::vector<Range<T>> res;

  auto start = problem.start();
  auto end = problem.end();

  if (end <= start) {
    throw XFDTDException("Invalid range");
  }

  auto size = end - start;

  for (int i = 0; i < num_procs; ++i) {
    auto task = decomposeSize(i, num_procs, size);
    res.emplace_back(task + start);
  }

  return res;
}

/**
 * @brief Decompose a task into sub-tasks. row major is the default.
 *
 * @param problem
 * @param type
 * @param nx
 * @param ny
 * @param nz
 * @param row_major
 * @return std::vector<Task<T>> A vector of sub-tasks
 */
template <typename T>
static auto decomposeTask(const Task<T>& problem, int nx, int ny, int nz,
                          bool row_major = true) -> std::vector<Task<T>> {
  if (nx <= 0 || ny <= 0 || nz <= 0) {
    return std::vector<Task<T>>{problem};
  }

  if (!row_major) {
    throw XFDTDException("Column major is not supported yet");
  }

  auto res = std::vector<Task<T>>{};
  auto x_ranges = decomposeRange(problem.xRange(), nx);
  auto y_ranges = decomposeRange(problem.yRange(), ny);
  auto z_ranges = decomposeRange(problem.zRange(), nz);
  for (auto i{0}; i < nx; ++i) {
    for (auto j{0}; j < ny; ++j) {
      for (auto k{0}; k < nz; ++k) {
        res.emplace_back(makeTask(x_ranges[i], y_ranges[j], z_ranges[k]));
      }
    }
  }

  return res;
}

}  // namespace xfdtd

#endif  // __XFDTD_CORE_DECOMPOSE_TASK_H__
