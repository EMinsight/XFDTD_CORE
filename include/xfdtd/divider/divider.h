#ifndef _XFDTD_CORE_DIVIDER_H_
#define _XFDTD_CORE_DIVIDER_H_

#include <xfdtd/exception/exception.h>

#include <cstddef>
#include <optional>
#include <sstream>
#include <vector>

namespace xfdtd {

class XFDTDDividerException : public XFDTDException {
 public:
  explicit XFDTDDividerException(
      const std::string& message = "XFDTD Divider Exception")
      : XFDTDException(message) {}
};

class Divider {
 public:
  enum class Type { UNDEFINED, X, Y, Z, XY, XZ, YZ, XYZ };

  template <typename T>
  class Range {
   public:
    Range() = default;
    Range(T start, T end) : _start(start), _end(end) {}

    auto operator[](std::size_t index) const {
      if (index == 0) {
        return _start;
      }
      return _end;
    }

    auto& operator[](std::size_t index) {
      if (index == 0) {
        return _start;
      }
      return _end;
    }

    auto operator+(T value) const {
      return Range{_start + value, _end + value};
    }

    auto operator-(T value) const {
      return Range{_start - value, _end - value};
    }

    auto start() const { return _start; }

    auto end() const { return _end; }

    auto size() const { return _end - _start; }

    auto valid() const { return _start < _end; }

    auto toString() const {
      std::stringstream ss;
      ss << "[" << _start << " " << _end << ")";
      return ss.str();
    }

   private:
    T _start{};
    T _end{};
  };

  template <typename T>
  struct Task {
    Range<T> _x_range;
    Range<T> _y_range;
    Range<T> _z_range;

    auto toString() const {
      std::stringstream ss;
      ss << "Task: ";
      ss << "x: " << _x_range.toString() << " ";
      ss << "y: " << _y_range.toString() << " ";
      ss << "z: " << _z_range.toString();
      return ss.str();
    }

    auto xRange() const { return _x_range; }

    auto yRange() const { return _y_range; }

    auto zRange() const { return _z_range; }

    auto valid() const {
      return _x_range.valid() && _y_range.valid() && _z_range.valid();
    }
  };

  using IndexRange = Range<std::size_t>;

  using IndexTask = Task<std::size_t>;

  template <typename T>
  static auto makeRange(T start, T end) {
    return Range<T>{start, end};
  }

  template <typename T>
  static auto makeTask(const Range<T>& x_range, const Range<T>& y_range,
                       const Range<T>& z_range) {
    return Task<T>{x_range, y_range, z_range};
  }

  static auto makeIndexRange(std::size_t start, std::size_t end) {
    return makeRange(start, end);
  }

  static auto makeIndexTask(const Range<std::size_t>& x_range,
                            const Range<std::size_t>& y_range,
                            const Range<std::size_t>& z_range) {
    return makeTask(x_range, y_range, z_range);
  }

  template <typename T>
  static bool intersected(const Range<T>& domain, const Range<T>& range) {
    return !(range[1] <= domain[0] || domain[1] <= range[0]);
  }

  template <typename T>
  static bool intersected(const Task<T>& task_1, const Task<T>& task_2) {
    return intersected(task_1._x_range, task_2._x_range) &&
           intersected(task_1._y_range, task_2._y_range) &&
           intersected(task_1._z_range, task_2._z_range);
  }

  template <typename T>
  static std::optional<Task<T>> taskIntersection(const Task<T>& task1,
                                                 const Task<T>& task2) {
    if (!intersected(task1, task2)) {
      return {};
    }

    auto x_range =
        Range<T>{std::max(task1.xRange().start(), task2.xRange().start()),
                 std::min(task1.xRange().end(), task2.xRange().end())};
    auto y_range =
        Range<T>{std::max(task1.yRange().start(), task2.yRange().start()),
                 std::min(task1.yRange().end(), task2.yRange().end())};
    auto z_range =
        Range<T>{std::max(task1.zRange().start(), task2.zRange().start()),
                 std::min(task1.zRange().end(), task2.zRange().end())};
    return Task<T>{x_range, y_range, z_range};
  }

  template <typename T>
  static auto divide(const Task<T>& problem, int num_procs, Type type) {
    if (num_procs <= 1) {
      return std::vector<Task<T>>{problem};
    }

    auto res = std::vector<Task<T>>{};
    switch (type) {
      case Type::X: {
        auto x_ranges = divide(problem.xRange(), num_procs);
        for (int i = 0; i < num_procs; ++i) {
          res.emplace_back(
              makeTask(x_ranges[i], problem.yRange(), problem.zRange()));
        }
        return res;
      }
      case Type::Y: {
        auto y_ranges = divide(problem.yRange(), num_procs);
        for (int i = 0; i < num_procs; ++i) {
          res.emplace_back(
              makeTask(problem.xRange(), y_ranges[i], problem.zRange()));
        }
        return res;
      }
      case Type::Z: {
        auto z_ranges = divide(problem.zRange(), num_procs);
        for (int i = 0; i < num_procs; ++i) {
          res.emplace_back(
              makeTask(problem.xRange(), problem.yRange(), z_ranges[i]));
        }
        return res;
      }
      case Type::XY: {
        if (num_procs == 2) {
          auto x_ranges = divide(problem.xRange(), 2);
          return std::vector<Task<T>>{
              makeTask(x_ranges[0], problem.yRange(), problem.zRange()),
              makeTask(x_ranges[1], problem.yRange(), problem.zRange())};
        }

        if (num_procs == 3) {
          auto x_ranges = divide(problem.xRange(), 2);
          auto l_y_ranges = divide(problem.yRange(), 2);
          return std::vector<Task<T>>{
              makeTask(x_ranges[0], l_y_ranges[0], problem.zRange()),
              makeTask(x_ranges[1], problem.yRange(), problem.zRange()),
              makeTask(x_ranges[0], l_y_ranges[1], problem.zRange())};
        }

        auto temp = divide<int>(makeRange<int>(0, num_procs), 4);
        auto x_ranges = divide(problem.xRange(), 2);
        auto y_ranges = divide(problem.yRange(), 2);
        for (int i = 0; i < 4; ++i) {
          auto x = i % 2;
          auto y = i / 2;
          auto unit =
              divide(makeIndexTask(x_ranges[x], y_ranges[y], problem.zRange()),
                     temp[i].size(), type);
          res.insert(res.end(), unit.begin(), unit.end());
        }

        return res;
      }
      default:
        throw XFDTDDividerException("Invalid type");
    }
  }

  template <typename T>
  static auto divide(const Range<T>& problem, int num_procs)
      -> std::vector<Range<T>> {
    std::vector<Range<T>> res;

    auto start = problem.start();
    auto end = problem.end();

    if (end <= start) {
      throw XFDTDDividerException("Invalid range");
    }

    auto size = end - start;

    for (int i = 0; i < num_procs; ++i) {
      auto task = divide(i, num_procs, size);
      res.push_back(task + start);
    }

    return res;
  }

  template <typename T>
  static Range<T> divide(int my_rank, int num_procs, const T& problem_size) {
    if (num_procs <= 0) {
      throw XFDTDDividerException(
          "Number of processes is less than or equal to zero");
    }

    if (problem_size < num_procs) {
      throw XFDTDDividerException(
          "Problem size is less than number of processes");
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
  auto fromTaskToType(const Task<T>& task) {
    throw XFDTDDividerException("Not implemented");
  }

 private:
};

}  // namespace xfdtd

#endif  // _XFDTD_CORE_DIVIDER_H_
